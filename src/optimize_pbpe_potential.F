!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_pbpe_potential
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE embed_types,                     ONLY: opt_pbpe_pot_type
   USE force_env_types,                 ONLY: force_env_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type,&
                                              set_mo_set
                                              
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_pbpe_potential'

   PUBLIC :: check_pbpe_system, prepare_pbpe_opt, release_pbpe_opt

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
! **************************************************************************************************
   SUBROUTINE check_pbpe_system(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      INTEGER                                            :: i_force_eval, n_force_eval, &
                                                            cluster_subsys_num, high_level_subsys_num, &
                                                            ref_subsys_num
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: naos, nmos
      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      cluster_subsys_num = 0
      high_level_subsys_num = 0
      ref_subsys_num = 0

      n_force_eval = SIZE(force_env%sub_force_env)
      IF (n_force_eval .NE. 4) THEN
         CPABORT("The number of systems is wrong.")
      END IF

      ALLOCATE (naos(n_force_eval))
      ALLOCATE (nmos(n_force_eval))

      ! Get info of subsystems 
      DO i_force_eval = 1, n_force_eval
         CALL get_qs_env(qs_env=force_env%sub_force_env(i_force_eval)%force_env%qs_env, &
                         do_kpoints=do_kpoints, dft_control=dft_control, mos=mos)
         CALL get_mo_set(mos(1), nao=naos(i_force_eval), nmo=nmos(i_force_eval))
         IF (dft_control%qs_control%cluster_embed_subsys) THEN
            cluster_subsys_num = i_force_eval
         ELSE IF (dft_control%qs_control%ref_embed_subsys) THEN
            ref_subsys_num = i_force_eval
         ELSE IF (dft_control%qs_control%high_level_embed_subsys) THEN
            high_level_subsys_num = i_force_eval
         END IF
         ! Check spin states
         IF (dft_control%nspins .NE. 1) THEN
            CPABORT("Open-shell PBPE not implemented.")
         END IF
         ! Check k-points
         IF ((dft_control%nimages .NE. 1) .OR. do_kpoints) THEN
            CPABORT("PBPE with k-points not implemented.")
         END IF
      END DO

      ! Check indices of different subsystems
      IF ((cluster_subsys_num .NE. 1) .AND. (cluster_subsys_num .NE. 2)) THEN
         CPABORT("Cluster subsystem must be the first or second FORCE_EVAL.")
      ELSE IF (ref_subsys_num .NE. 3) THEN
         CPABORT("Reference subsystem must be the third FORCE_EVAL.")
      ELSE IF (high_level_subsys_num .NE. 4) THEN
         CPABORT("High-level subsystem must be the fourth FORCE_EVAL.")
      END IF 

      ! Check numbers of AOs in different subsystems
      !IF (naos(1) + naos(2) .NE. naos(3)) THEN
      !   CPABORT("The sum of nAOs in two low-level subsystem does not equal to nAOs in the total system.")
      !ELSE IF (naos(2) .NE. naos(4)) THEN
      !   CPABORT("The low-level and high-level calculations should use the same basis set.")
      !END IF

      DEALLOCATE (naos)
      DEALLOCATE (nmos)

   END SUBROUTINE check_pbpe_system

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE init_pbpe_matrices(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL get_qs_env(qs_env=force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env, &
                      para_env=para_env)

      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      NULLIFY (opt_pbpe%dm_1, opt_pbpe%dm_2, opt_pbpe%dm_total)
      NULLIFY (opt_pbpe%dm_1_beta, opt_pbpe%dm_2_beta, opt_pbpe%dm_total_beta)
      NULLIFY (opt_pbpe%dm_1_diff, opt_pbpe%dm_2_diff, opt_pbpe%dm_total_diff)
      NULLIFY (opt_pbpe%dm_1_diff_beta, opt_pbpe%dm_2_diff_beta, opt_pbpe%dm_total_diff_beta)
      NULLIFY (opt_pbpe%dm_cluster, opt_pbpe%dm_environment)
      NULLIFY (opt_pbpe%dm_cluster_beta, opt_pbpe%dm_environment_beta)
      NULLIFY (opt_pbpe%dm_cluster_diff, opt_pbpe%dm_environment_diff)
      NULLIFY (opt_pbpe%dm_cluster_diff_beta, opt_pbpe%dm_environment_diff_beta)
      NULLIFY (opt_pbpe%S_total, opt_pbpe%KS_total)
      NULLIFY (opt_pbpe%S_1, opt_pbpe%S_2)
      NULLIFY (opt_pbpe%KS_1, opt_pbpe%KS_2)
      NULLIFY (opt_pbpe%KS_1_in_2, opt_pbpe%KS_2_in_1)
      NULLIFY (opt_pbpe%S_12, opt_pbpe%S_21)
      NULLIFY (opt_pbpe%KS_12, opt_pbpe%KS_21)
      NULLIFY (opt_pbpe%work1_1, opt_pbpe%work1_2)
      NULLIFY (opt_pbpe%work2_1, opt_pbpe%work2_2)
      NULLIFY (opt_pbpe%ortho_1, opt_pbpe%ortho_2)
      NULLIFY (opt_pbpe%P_1, opt_pbpe%P_2)
      NULLIFY (opt_pbpe%P_environment) 

      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(3), ncol_global=opt_pbpe%naos_in_pbpe(3))
      ALLOCATE (opt_pbpe%dm_total, opt_pbpe%dm_total_diff, opt_pbpe%S_total, opt_pbpe%KS_total)
      CALL cp_fm_create(opt_pbpe%dm_total, fm_struct, name="dm_total")
      CALL cp_fm_create(opt_pbpe%dm_total_diff, fm_struct, name="dm_total_diff")
      CALL cp_fm_create(opt_pbpe%S_total, fm_struct, name="S_total")
      CALL cp_fm_create(opt_pbpe%KS_total, fm_struct, name="KS_total")

      CALL cp_fm_set_all(opt_pbpe%dm_total, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_total_diff, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%S_total, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_total, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_pbpe%open_shell_embed) THEN
         ALLOCATE (opt_pbpe%dm_total_beta, opt_pbpe%dm_total_diff_beta)
         CALL cp_fm_create(opt_pbpe%dm_total_beta, fm_struct, name="dm_total_beta")
         CALL cp_fm_create(opt_pbpe%dm_total_diff_beta, fm_struct, name="dm_total_diff_beta")

         CALL cp_fm_set_all(opt_pbpe%dm_total_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_total_diff_beta, 0.0_dp)
      END IF

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(1), ncol_global=opt_pbpe%naos_in_pbpe(1))
      ALLOCATE (opt_pbpe%dm_1, opt_pbpe%dm_1_diff, opt_pbpe%S_1, opt_pbpe%KS_1, opt_pbpe%KS_1_in_2, & 
                opt_pbpe%work1_1, opt_pbpe%work2_1, opt_pbpe%ortho_1, opt_pbpe%P_2)
      CALL cp_fm_create(opt_pbpe%dm_1, fm_struct, name="dm_1")
      CALL cp_fm_create(opt_pbpe%dm_1_diff, fm_struct, name="dm_1_diff")
      CALL cp_fm_create(opt_pbpe%S_1, fm_struct, name="S_1")
      CALL cp_fm_create(opt_pbpe%KS_1, fm_struct, name="KS_1")
      CALL cp_fm_create(opt_pbpe%KS_1_in_2, fm_struct, name="KS_1_in_2")
      CALL cp_fm_create(opt_pbpe%work1_1, fm_struct, name="work1_1")
      CALL cp_fm_create(opt_pbpe%work2_1, fm_struct, name="work2_1")
      CALL cp_fm_create(opt_pbpe%ortho_1, fm_struct, name="ortho_1")
      CALL cp_fm_create(opt_pbpe%P_2, fm_struct, name="P_2")

      CALL cp_fm_set_all(opt_pbpe%dm_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_1_diff, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%S_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_1_in_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work1_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work2_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ortho_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%P_2, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_pbpe%open_shell_embed) THEN
         ALLOCATE (opt_pbpe%dm_1_beta, opt_pbpe%dm_1_diff_beta)
         CALL cp_fm_create(opt_pbpe%dm_1_beta, fm_struct, name="dm_1_beta")
         CALL cp_fm_create(opt_pbpe%dm_1_diff_beta, fm_struct, name="dm_1_diff_beta")

         CALL cp_fm_set_all(opt_pbpe%dm_1_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_1_diff_beta, 0.0_dp)
      END IF

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(2), ncol_global=opt_pbpe%naos_in_pbpe(2))
      ALLOCATE (opt_pbpe%dm_2, opt_pbpe%dm_2_diff, opt_pbpe%S_2, opt_pbpe%KS_2, opt_pbpe%KS_2_in_1, &
                opt_pbpe%work1_2, opt_pbpe%work2_2, opt_pbpe%ortho_2, opt_pbpe%P_1)
      CALL cp_fm_create(opt_pbpe%dm_2, fm_struct, name="dm_2")
      CALL cp_fm_create(opt_pbpe%dm_2_diff, fm_struct, name="dm_2_diff")
      CALL cp_fm_create(opt_pbpe%S_2, fm_struct, name="S_2")
      CALL cp_fm_create(opt_pbpe%KS_2, fm_struct, name="KS_2")
      CALL cp_fm_create(opt_pbpe%KS_2_in_1, fm_struct, name="KS_2_in_1")
      CALL cp_fm_create(opt_pbpe%work1_2, fm_struct, name="work1_2")
      CALL cp_fm_create(opt_pbpe%work2_2, fm_struct, name="work2_2")
      CALL cp_fm_create(opt_pbpe%ortho_2, fm_struct, name="ortho_2")
      CALL cp_fm_create(opt_pbpe%P_1, fm_struct, name="P_1")

      CALL cp_fm_set_all(opt_pbpe%dm_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_2_diff, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%S_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_2_in_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work1_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work2_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ortho_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%P_1, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_pbpe%open_shell_embed) THEN
         ALLOCATE (opt_pbpe%dm_2_beta, opt_pbpe%dm_2_diff_beta)
         CALL cp_fm_create(opt_pbpe%dm_2_beta, fm_struct, name="dm_2_beta")
         CALL cp_fm_create(opt_pbpe%dm_2_diff_beta, fm_struct, name="dm_2_diff_beta")

         CALL cp_fm_set_all(opt_pbpe%dm_2_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_2_diff_beta, 0.0_dp)
      END IF

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(1), ncol_global=opt_pbpe%naos_in_pbpe(2))
      ALLOCATE (opt_pbpe%S_12, opt_pbpe%KS_12)
      CALL cp_fm_create(opt_pbpe%S_12, fm_struct, name="S_12")
      CALL cp_fm_create(opt_pbpe%KS_12, fm_struct, name="KS_12")

      CALL cp_fm_set_all(opt_pbpe%S_12, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_12, 0.0_dp)

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(2), ncol_global=opt_pbpe%naos_in_pbpe(1))
      ALLOCATE (opt_pbpe%S_21, opt_pbpe%KS_21)
      CALL cp_fm_create(opt_pbpe%S_21, fm_struct, name="S_21")
      CALL cp_fm_create(opt_pbpe%KS_21, fm_struct, name="KS_21")

      CALL cp_fm_set_all(opt_pbpe%S_21, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%KS_21, 0.0_dp)

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      ! Release context
      CALL cp_blacs_env_release(blacs_env)



      ! Associate matrix pointers
      IF (opt_pbpe%cluster_subsys_num .EQ. 1) THEN
         opt_pbpe%dm_cluster => opt_pbpe%dm_1
         opt_pbpe%dm_environment => opt_pbpe%dm_2
         opt_pbpe%dm_cluster_diff => opt_pbpe%dm_1_diff
         opt_pbpe%dm_environment_diff => opt_pbpe%dm_2_diff
         ! Beta spin counterparts
         IF (opt_pbpe%open_shell_embed) THEN
            opt_pbpe%dm_cluster_beta => opt_pbpe%dm_1_beta
            opt_pbpe%dm_environment_beta => opt_pbpe%dm_2_beta
            opt_pbpe%dm_cluster_diff_beta => opt_pbpe%dm_1_diff_beta
            opt_pbpe%dm_environment_diff_beta => opt_pbpe%dm_2_diff_beta
         END IF
      ELSE IF (opt_pbpe%cluster_subsys_num .EQ. 2) THEN
         opt_pbpe%dm_cluster => opt_pbpe%dm_2
         opt_pbpe%dm_environment => opt_pbpe%dm_1
         opt_pbpe%dm_cluster_diff => opt_pbpe%dm_2_diff
         opt_pbpe%dm_environment_diff => opt_pbpe%dm_1_diff
         ! Beta spin counterparts
         IF (opt_pbpe%open_shell_embed) THEN
            opt_pbpe%dm_cluster_beta => opt_pbpe%dm_2_beta
            opt_pbpe%dm_environment_beta => opt_pbpe%dm_1_beta
            opt_pbpe%dm_cluster_diff_beta => opt_pbpe%dm_2_diff_beta
            opt_pbpe%dm_environment_diff_beta => opt_pbpe%dm_1_diff_beta
         END IF
      END IF

   END SUBROUTINE init_pbpe_matrices

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE init_pbpe_mo_sets(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      INTEGER                                            :: i_force_eval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_1, mo_coeff_2
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL get_qs_env(qs_env=force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env, &
                      para_env=para_env)

      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      NULLIFY (opt_pbpe%mos_1, opt_pbpe%mos_2, opt_pbpe%mos_cluster, opt_pbpe%mos_environment)

      ALLOCATE (opt_pbpe%mos_1, opt_pbpe%mos_2)

      NULLIFY (mo_coeff_1, mo_coeff_2)
      NULLIFY (fm_struct)

      ! Initialize mos_1
      CALL get_qs_env(qs_env=force_env%sub_force_env(1)%force_env%qs_env, mos=mos)
      CALL duplicate_mo_set(opt_pbpe%mos_1, mos(1))
      ! Copy mos_1%mo_coeff to mo_coeff_1
      CALL get_mo_set(opt_pbpe%mos_1, mo_coeff=mo_coeff)
      ALLOCATE (mo_coeff_1)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(1), ncol_global=opt_pbpe%nmos(1))
      CALL cp_fm_create(mo_coeff_1, fm_struct)
      CALL cp_fm_set_all(mo_coeff_1, 0.0_dp)
      IF (.NOT. opt_pbpe%abs_localization) THEN
         CALL cp_fm_to_fm_submat(mo_coeff, mo_coeff_1, opt_pbpe%naos(1), opt_pbpe%nmos(1), 1, 1, 1, 1)
      ELSE
         CALL cp_fm_to_fm(mo_coeff, mo_coeff_1)
      END IF
      ! Replace mos_1%mo_coeff with mo_coeff_1
      CALL set_mo_set(opt_pbpe%mos_1, nao=opt_pbpe%naos_in_pbpe(1))
      CALL cp_fm_release(opt_pbpe%mos_1%mo_coeff)
      DEALLOCATE (opt_pbpe%mos_1%mo_coeff)
      opt_pbpe%mos_1%mo_coeff => mo_coeff_1

      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      ! Initialize mos_2
      CALL get_qs_env(qs_env=force_env%sub_force_env(2)%force_env%qs_env, mos=mos)
      CALL duplicate_mo_set(opt_pbpe%mos_2, mos(1))
      ! Copy mos_2%mo_coeff to mo_coeff_2
      CALL get_mo_set(opt_pbpe%mos_2, mo_coeff=mo_coeff)
      ALLOCATE (mo_coeff_2)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(2), ncol_global=opt_pbpe%nmos(2))
      CALL cp_fm_create(mo_coeff_2, fm_struct)
      CALL cp_fm_set_all(mo_coeff_2, 0.0_dp)
      IF (.NOT. opt_pbpe%abs_localization) THEN
         CALL cp_fm_to_fm_submat(mo_coeff, mo_coeff_2, opt_pbpe%naos(2), opt_pbpe%nmos(2), 1, 1, opt_pbpe%naos(1) + 1, 1)
      ELSE
         CALL cp_fm_to_fm(mo_coeff, mo_coeff_2)
      END IF
      ! Replace mos_2%mo_coeff with mo_coeff_2
      CALL set_mo_set(opt_pbpe%mos_2, nao=opt_pbpe%naos_in_pbpe(2))
      CALL cp_fm_release(opt_pbpe%mos_2%mo_coeff)
      DEALLOCATE (opt_pbpe%mos_2%mo_coeff)
      opt_pbpe%mos_2%mo_coeff => mo_coeff_2

      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_blacs_env_release(blacs_env)

      ! Associate mos pointers
      IF (opt_pbpe%cluster_subsys_num .EQ. 1) THEN
         opt_pbpe%mos_cluster => opt_pbpe%mos_1
         opt_pbpe%mos_environment => opt_pbpe%mos_2
      ELSE IF (opt_pbpe%cluster_subsys_num .EQ. 2) THEN
         opt_pbpe%mos_cluster => opt_pbpe%mos_2
         opt_pbpe%mos_environment => opt_pbpe%mos_1
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
!> \param opt_pbpe_section ...
! **************************************************************************************************
   SUBROUTINE prepare_pbpe_opt(force_env, opt_pbpe, opt_pbpe_section)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      TYPE(section_vals_type), POINTER                   :: opt_pbpe_section

      INTEGER                                            :: i_force_eval, n_force_eval
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      ! Read the input
      CALL read_opt_pbpe_section(opt_pbpe, opt_pbpe_section)

      n_force_eval = SIZE(force_env%sub_force_env)

      ALLOCATE (opt_pbpe%naos(n_force_eval))
      ALLOCATE (opt_pbpe%naos_in_pbpe(n_force_eval))
      ALLOCATE (opt_pbpe%nmos(n_force_eval))

      ! Get info of subsystems 
      DO i_force_eval = 1, n_force_eval
         CALL get_qs_env(qs_env=force_env%sub_force_env(i_force_eval)%force_env%qs_env, &
                         dft_control=dft_control, mos=mos)
         CALL get_mo_set(mos(1), nao=opt_pbpe%naos(i_force_eval), nmo=opt_pbpe%nmos(i_force_eval))
         IF (dft_control%qs_control%cluster_embed_subsys) THEN
            opt_pbpe%cluster_subsys_num = i_force_eval
         ELSE IF (dft_control%qs_control%ref_embed_subsys) THEN
            opt_pbpe%ref_subsys_num = i_force_eval
         ELSE IF (dft_control%qs_control%high_level_embed_subsys) THEN
            opt_pbpe%high_level_subsys_num = i_force_eval
         END IF
      END DO

      IF (opt_pbpe%cluster_subsys_num .EQ. 1) THEN
         opt_pbpe%environment_subsys_num = 2
      ELSE IF (opt_pbpe%cluster_subsys_num .EQ. 2) THEN
         opt_pbpe%environment_subsys_num = 1
      END IF

      DO i_force_eval = 1, n_force_eval
         opt_pbpe%naos_in_pbpe(i_force_eval) = opt_pbpe%naos(i_force_eval) 
      END DO

      ! modify naos_in_pbpe if absolute localization is not used
      IF (.NOT. opt_pbpe%abs_localization) THEN
         opt_pbpe%naos_in_pbpe(opt_pbpe%cluster_subsys_num) = opt_pbpe%naos_in_pbpe(opt_pbpe%ref_subsys_num)
         opt_pbpe%naos_in_pbpe(opt_pbpe%environment_subsys_num) = opt_pbpe%naos_in_pbpe(opt_pbpe%ref_subsys_num)
         opt_pbpe%naos_in_pbpe(opt_pbpe%high_level_subsys_num) = opt_pbpe%naos_in_pbpe(opt_pbpe%ref_subsys_num)
      END IF

      ! Initialize MO sets
      CALL init_pbpe_mo_sets(force_env, opt_pbpe)

      ! Initialize cp_fm matrices
      CALL init_pbpe_matrices(force_env, opt_pbpe)

   END SUBROUTINE prepare_pbpe_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE read_opt_pbpe_section(opt_pbpe, opt_pbpe_section)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      TYPE(section_vals_type), POINTER                   :: opt_pbpe_section

      ! Read keywords

      CALL section_vals_val_get(opt_pbpe_section, "N_ITER", &
                                i_val=opt_pbpe%n_iter)

      CALL section_vals_val_get(opt_pbpe_section, "ABS_LOC", &
                                l_val=opt_pbpe%abs_localization)

      CALL section_vals_val_get(opt_pbpe_section, "DM_CONV_MAX", &
                                r_val=opt_pbpe%conv_max)

      CALL section_vals_val_get(opt_pbpe_section, "DM_CONV_RMSE", &
                                r_val=opt_pbpe%conv_rmse)

   END SUBROUTINE read_opt_pbpe_section

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE release_pbpe_matrices(opt_pbpe)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      NULLIFY (opt_pbpe%dm_cluster)
      NULLIFY (opt_pbpe%dm_environment)
      NULLIFY (opt_pbpe%dm_cluster_diff)
      NULLIFY (opt_pbpe%dm_environment_diff)

      IF (opt_pbpe%open_shell_embed) THEN
         NULLIFY (opt_pbpe%dm_cluster_beta)
         NULLIFY (opt_pbpe%dm_environment_beta)
         NULLIFY (opt_pbpe%dm_cluster_diff_beta)
         NULLIFY (opt_pbpe%dm_environment_diff_beta)
      END IF

      IF (ASSOCIATED(opt_pbpe%dm_total)) THEN
         CALL cp_fm_release(opt_pbpe%dm_total)
         DEALLOCATE (opt_pbpe%dm_total)
         NULLIFY (opt_pbpe%dm_total)
      END IF
      IF (ASSOCIATED(opt_pbpe%dm_total_diff)) THEN
         CALL cp_fm_release(opt_pbpe%dm_total_diff)
         DEALLOCATE (opt_pbpe%dm_total_diff)
         NULLIFY (opt_pbpe%dm_total_diff)
      END IF
      IF (ASSOCIATED(opt_pbpe%S_total)) THEN
         CALL cp_fm_release(opt_pbpe%S_total)
         DEALLOCATE (opt_pbpe%S_total)
         NULLIFY (opt_pbpe%S_total)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_total)) THEN
         CALL cp_fm_release(opt_pbpe%KS_total)
         DEALLOCATE (opt_pbpe%KS_total)
         NULLIFY (opt_pbpe%KS_total)
      END IF 

      IF (opt_pbpe%open_shell_embed) THEN
         IF (ASSOCIATED(opt_pbpe%dm_total_beta)) THEN
            CALL cp_fm_release(opt_pbpe%dm_total_beta)
            DEALLOCATE (opt_pbpe%dm_total_beta)
            NULLIFY (opt_pbpe%dm_total_beta)
         END IF
         IF (ASSOCIATED(opt_pbpe%dm_total_diff_beta)) THEN
            CALL cp_fm_release(opt_pbpe%dm_total_diff_beta)
            DEALLOCATE (opt_pbpe%dm_total_diff_beta)
            NULLIFY (opt_pbpe%dm_total_diff_beta)
         END IF
      END IF

      IF (ASSOCIATED(opt_pbpe%dm_1)) THEN
         CALL cp_fm_release(opt_pbpe%dm_1)
         DEALLOCATE (opt_pbpe%dm_1)
         NULLIFY (opt_pbpe%dm_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%dm_1_diff)) THEN
         CALL cp_fm_release(opt_pbpe%dm_1_diff)
         DEALLOCATE (opt_pbpe%dm_1_diff)
         NULLIFY (opt_pbpe%dm_1_diff)
      END IF
      IF (ASSOCIATED(opt_pbpe%S_1)) THEN
         CALL cp_fm_release(opt_pbpe%S_1)
         DEALLOCATE (opt_pbpe%S_1)
         NULLIFY (opt_pbpe%S_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_1)) THEN
         CALL cp_fm_release(opt_pbpe%KS_1)
         DEALLOCATE (opt_pbpe%KS_1)
         NULLIFY (opt_pbpe%KS_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_1_in_2)) THEN
         CALL cp_fm_release(opt_pbpe%KS_1_in_2)
         DEALLOCATE (opt_pbpe%KS_1_in_2)
         NULLIFY (opt_pbpe%KS_1_in_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%work1_1)) THEN
         CALL cp_fm_release(opt_pbpe%work1_1)
         DEALLOCATE (opt_pbpe%work1_1)
         NULLIFY (opt_pbpe%work1_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%work2_1)) THEN
         CALL cp_fm_release(opt_pbpe%work2_1)
         DEALLOCATE (opt_pbpe%work2_1)
         NULLIFY (opt_pbpe%work2_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%ortho_1)) THEN
         CALL cp_fm_release(opt_pbpe%ortho_1)
         DEALLOCATE (opt_pbpe%ortho_1)
         NULLIFY (opt_pbpe%ortho_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%P_2)) THEN
         CALL cp_fm_release(opt_pbpe%P_2)
         DEALLOCATE (opt_pbpe%P_2)
         NULLIFY (opt_pbpe%P_2)
      END IF

      IF (opt_pbpe%open_shell_embed) THEN
         IF (ASSOCIATED(opt_pbpe%dm_1_beta)) THEN
            CALL cp_fm_release(opt_pbpe%dm_1_beta)
            DEALLOCATE (opt_pbpe%dm_1_beta)
            NULLIFY (opt_pbpe%dm_1_beta)
         END IF
         IF (ASSOCIATED(opt_pbpe%dm_1_diff_beta)) THEN
            CALL cp_fm_release(opt_pbpe%dm_1_diff_beta)
            DEALLOCATE (opt_pbpe%dm_1_diff_beta)
            NULLIFY (opt_pbpe%dm_1_diff_beta)
         END IF
      END IF

      IF (ASSOCIATED(opt_pbpe%dm_2)) THEN
         CALL cp_fm_release(opt_pbpe%dm_2)
         DEALLOCATE (opt_pbpe%dm_2)
         NULLIFY (opt_pbpe%dm_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%dm_2_diff)) THEN
         CALL cp_fm_release(opt_pbpe%dm_2_diff)
         DEALLOCATE (opt_pbpe%dm_2_diff)
         NULLIFY (opt_pbpe%dm_2_diff)
      END IF
      IF (ASSOCIATED(opt_pbpe%S_2)) THEN
         CALL cp_fm_release(opt_pbpe%S_2)
         DEALLOCATE (opt_pbpe%S_2)
         NULLIFY (opt_pbpe%S_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_2)) THEN
         CALL cp_fm_release(opt_pbpe%KS_2)
         DEALLOCATE (opt_pbpe%KS_2)
         NULLIFY (opt_pbpe%KS_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_2_in_1)) THEN
         CALL cp_fm_release(opt_pbpe%KS_2_in_1)
         DEALLOCATE (opt_pbpe%KS_2_in_1)
         NULLIFY (opt_pbpe%KS_2_in_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%work1_2)) THEN
         CALL cp_fm_release(opt_pbpe%work1_2)
         DEALLOCATE (opt_pbpe%work1_2)
         NULLIFY (opt_pbpe%work1_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%work2_2)) THEN
         CALL cp_fm_release(opt_pbpe%work2_2)
         DEALLOCATE (opt_pbpe%work2_2)
         NULLIFY (opt_pbpe%work2_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%ortho_2)) THEN
         CALL cp_fm_release(opt_pbpe%ortho_2)
         DEALLOCATE (opt_pbpe%ortho_2)
         NULLIFY (opt_pbpe%ortho_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%P_1)) THEN
         CALL cp_fm_release(opt_pbpe%P_1)
         DEALLOCATE (opt_pbpe%P_1)
         NULLIFY (opt_pbpe%P_1)
      END IF

      IF (opt_pbpe%open_shell_embed) THEN
         IF (ASSOCIATED(opt_pbpe%dm_2_beta)) THEN
            CALL cp_fm_release(opt_pbpe%dm_2_beta)
            DEALLOCATE (opt_pbpe%dm_2_beta)
            NULLIFY (opt_pbpe%dm_2_beta)
         END IF
         IF (ASSOCIATED(opt_pbpe%dm_2_diff_beta)) THEN
            CALL cp_fm_release(opt_pbpe%dm_2_diff_beta)
            DEALLOCATE (opt_pbpe%dm_2_diff_beta)
            NULLIFY (opt_pbpe%dm_2_diff_beta)
         END IF
      END IF

      IF (ASSOCIATED(opt_pbpe%S_12)) THEN
         CALL cp_fm_release(opt_pbpe%S_12)
         DEALLOCATE (opt_pbpe%S_12)
         NULLIFY (opt_pbpe%S_12)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_12)) THEN
         CALL cp_fm_release(opt_pbpe%KS_12)
         DEALLOCATE (opt_pbpe%KS_12)
         NULLIFY (opt_pbpe%KS_12)
      END IF

      IF (ASSOCIATED(opt_pbpe%S_21)) THEN
         CALL cp_fm_release(opt_pbpe%S_21)
         DEALLOCATE (opt_pbpe%S_21)
         NULLIFY (opt_pbpe%S_21)
      END IF
      IF (ASSOCIATED(opt_pbpe%KS_21)) THEN
         CALL cp_fm_release(opt_pbpe%KS_21)
         DEALLOCATE (opt_pbpe%KS_21)
         NULLIFY (opt_pbpe%KS_21)
      END IF

   END SUBROUTINE release_pbpe_matrices

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE release_pbpe_mo_sets(opt_pbpe)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      NULLIFY (opt_pbpe%mos_cluster)
      NULLIFY (opt_pbpe%mos_environment)

      IF (ASSOCIATED(opt_pbpe%mos_1)) THEN
         CALL deallocate_mo_set(opt_pbpe%mos_1)
         DEALLOCATE (opt_pbpe%mos_1)
         NULLIFY (opt_pbpe%mos_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%mos_2)) THEN
         CALL deallocate_mo_set(opt_pbpe%mos_2)
         DEALLOCATE (opt_pbpe%mos_2)
         NULLIFY (opt_pbpe%mos_2)
      END IF

   END SUBROUTINE release_pbpe_mo_sets

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE release_pbpe_opt(opt_pbpe)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CALL release_pbpe_matrices(opt_pbpe)

      CALL release_pbpe_mo_sets(opt_pbpe)

      DEALLOCATE (opt_pbpe%naos)
      DEALLOCATE (opt_pbpe%naos_in_pbpe)
      DEALLOCATE (opt_pbpe%nmos)

   END SUBROUTINE release_pbpe_opt

END MODULE optimize_pbpe_potential