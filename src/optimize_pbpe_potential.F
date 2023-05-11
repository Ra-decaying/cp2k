!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_pbpe_potential
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_dbcsr_output,                 ONLY: write_fm_with_basis_info
   USE cp_fm_basic_linalg,              ONLY: cp_fm_geadd,&
                                              cp_fm_norm,&
                                              cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE embed_types,                     ONLY: opt_pbpe_pot_type
   USE force_env_types,                 ONLY: force_env_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type,&
                                              set_mo_set
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_loop_utils,               ONLY: qs_scf_rho_update
   USE qs_scf_methods,                  ONLY: eigensolver_symm
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_pbpe_potential'

                                              
   PUBLIC :: build_pbpe_subsys_dm, build_pbpe_ref_subsys_dm, check_pbpe_system, &
             calculate_embedded_ks_matrix, calculate_pbpe_dm_difference, &
             check_pbpe_low_level_convergence, copy_ks_to_pbpe, &
             pbpe_do_diagonalization, pbpe_initialization_print, &
             pbpe_low_level_loop_print, pbpe_summary_print, prepare_pbpe_opt, &
             release_pbpe_opt, update_pbpe_density_and_ks_qs_env

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
!> \param subsys_num ...
!> \param density_matrix ...
! **************************************************************************************************
   SUBROUTINE build_pbpe_subsys_dm(opt_pbpe, subsys_num, density_matrix)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      INTEGER, INTENT(IN)                                :: subsys_num
      TYPE(cp_fm_type), POINTER                          :: density_matrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_pbpe_subsys_dm'

      INTEGER                                            :: homo, nao
      REAL(KIND=dp)                                      :: coeff
      TYPE(cp_fm_type), POINTER                          :: mo_coeff

      CPASSERT((subsys_num .EQ. 1) .OR. (subsys_num .EQ. 2))

      coeff = 2.0_dp
      IF (subsys_num .EQ. 1) THEN
         CALL get_mo_set(mo_set=opt_pbpe%mos_1(1), mo_coeff=mo_coeff, nao=nao, homo=homo)
         CALL parallel_gemm(transa="N", transb="T", m=nao, n=nao, k=homo, &
                            alpha=coeff, matrix_a=mo_coeff, matrix_b=mo_coeff, &
                            beta=0.0_dp, matrix_c=density_matrix)
      ELSE IF (subsys_num .EQ. 2) THEN
         CALL get_mo_set(mo_set=opt_pbpe%mos_2(1), mo_coeff=mo_coeff, nao=nao, homo=homo)
         CALL parallel_gemm(transa="N", transb="T", m=nao, n=nao, k=homo, &
                            alpha=coeff, matrix_a=mo_coeff, matrix_b=mo_coeff, &
                            beta=0.0_dp, matrix_c=density_matrix)
      END IF

   END SUBROUTINE build_pbpe_subsys_dm

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
!> \param density_matrix ...
!> \param use_new_subsys_dm ...
! **************************************************************************************************
   SUBROUTINE build_pbpe_ref_subsys_dm(opt_pbpe, density_matrix, use_new_subsys_dm)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      TYPE(cp_fm_type), POINTER                          :: density_matrix
      LOGICAL                                            :: use_new_subsys_dm

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_pbpe_ref_sys_dm'

      TYPE(cp_fm_type), POINTER                          :: dm_1, dm_2

      IF (use_new_subsys_dm) THEN
         dm_1 => opt_pbpe%dm_1_new
         dm_2 => opt_pbpe%dm_2_new
      ELSE
         dm_1 => opt_pbpe%dm_1
         dm_2 => opt_pbpe%dm_2
      ENDIF

      IF (opt_pbpe%abs_localization) THEN
         CALL cp_fm_to_fm_submat(dm_1, density_matrix, opt_pbpe%naos(1), opt_pbpe%naos(1), &
                                 1, 1, 1, 1)
         CALL cp_fm_to_fm_submat(dm_2, density_matrix, opt_pbpe%naos(2), opt_pbpe%naos(2), &
                                 1, 1, opt_pbpe%naos(2) + 1, opt_pbpe%naos(2) + 1)
      ELSE
         CALL cp_fm_to_fm(opt_pbpe%dm_1, density_matrix)
         CALL cp_fm_scale_and_add(1.0_dp, density_matrix, 1.0_dp, opt_pbpe%dm_2)
      END IF

   END SUBROUTINE build_pbpe_ref_subsys_dm

! **************************************************************************************************
!> \brief
!> \param opt_pbpe ...
!> \param subsys_num ...
! **************************************************************************************************
   SUBROUTINE calculate_embedded_ks_matrix(opt_pbpe, subsys_num)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      INTEGER, INTENT(IN)                                :: subsys_num

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calculate_embedded_ks_matrix'

      CPASSERT((subsys_num .EQ. 1) .OR. (subsys_num .EQ. 2))
      IF (subsys_num .EQ. 1) THEN
         ! fm_tmp1_1 = ks_12 - fermi_energy * s_12
         CALL cp_fm_to_fm(opt_pbpe%ks_12, opt_pbpe%fm_tmp1_1)
         CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%fm_tmp1_1, -1.0_dp * opt_pbpe%fermi_energy, opt_pbpe%s_12)
         ! fm_tmp2_1 = fm_tmp1_1 * dm_2
         CALL parallel_gemm('N', 'N', opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(2), &
                            1.0_dp, opt_pbpe%fm_tmp1_1, opt_pbpe%dm_2, 0.0_dp, opt_pbpe%fm_tmp2_1)
         ! fm_tmp3_1 = fm_tmp2_1 * s_21
         CALL parallel_gemm('N', 'N', opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(2), &
                            1.0_dp, opt_pbpe%fm_tmp2_1, opt_pbpe%s_21, 0.0_dp, opt_pbpe%fm_tmp3_1)
         ! projection_2 = -0.5 * fm_tmp3_1 * (fm_tmp3_1 ^ T)
         CALL cp_fm_to_fm(opt_pbpe%fm_tmp3_1, opt_pbpe%projection_2)
         CALL cp_fm_geadd(-0.5_dp, 'T', opt_pbpe%fm_tmp3_1, -0.5_dp, opt_pbpe%projection_2)
         ! Calculate ks_1_in_2
         CALL cp_fm_to_fm(opt_pbpe%ks_1, opt_pbpe%ks_1_in_2)
         CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%ks_1_in_2, 1.0_dp, opt_pbpe%projection_2)
      ELSE IF (subsys_num .EQ. 2) THEN
         ! fm_tmp1_2 = ks_21 - fermi_energy * s_21
         CALL cp_fm_to_fm(opt_pbpe%ks_21, opt_pbpe%fm_tmp1_2)
         CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%fm_tmp1_2, -1.0_dp * opt_pbpe%fermi_energy, opt_pbpe%s_21)
         ! fm_tmp2_2 = fm_tmp1_2 * dm_1
         CALL parallel_gemm('N', 'N', opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(1), &
                            1.0_dp, opt_pbpe%fm_tmp1_2, opt_pbpe%dm_1, 0.0_dp, opt_pbpe%fm_tmp2_2)
         ! fm_tmp3_2 = fm_tmp2_2 * s_12
         CALL parallel_gemm('N', 'N', opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(1), &
                            1.0_dp, opt_pbpe%fm_tmp2_2, opt_pbpe%s_12, 0.0_dp, opt_pbpe%fm_tmp3_2)
         ! projection_1 = -0.5 * fm_tmp3_2 * (fm_tmp3_2 ^ T)
         CALL cp_fm_to_fm(opt_pbpe%fm_tmp3_2, opt_pbpe%projection_1)
         CALL cp_fm_geadd(-0.5_dp, 'T', opt_pbpe%fm_tmp3_2, -0.5_dp, opt_pbpe%projection_1)
         ! Calculate ks_2_in_1
         CALL cp_fm_to_fm(opt_pbpe%ks_2, opt_pbpe%ks_2_in_1)
         CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%ks_2_in_1, 1.0_dp, opt_pbpe%projection_1)
      END IF

   END SUBROUTINE calculate_embedded_ks_matrix

! **************************************************************************************************
!> \brief
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE calculate_pbpe_dm_difference(opt_pbpe)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calculate_pbpe_dm_difference'

      CALL cp_fm_to_fm(opt_pbpe%dm_1_new, opt_pbpe%dm_1_diff)
      CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%dm_1_diff, -1.0_dp, opt_pbpe%dm_1)
      CALL cp_fm_to_fm(opt_pbpe%dm_2_new, opt_pbpe%dm_2_diff)
      CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%dm_2_diff, -1.0_dp, opt_pbpe%dm_2)
      CALL cp_fm_to_fm(opt_pbpe%dm_total_new, opt_pbpe%dm_total_diff)
      CALL cp_fm_scale_and_add(1.0_dp, opt_pbpe%dm_total_diff, -1.0_dp, opt_pbpe%dm_total)

      CALL cp_fm_to_fm(opt_pbpe%dm_1_new, opt_pbpe%dm_1)
      CALL cp_fm_to_fm(opt_pbpe%dm_2_new, opt_pbpe%dm_2)
      CALL cp_fm_to_fm(opt_pbpe%dm_total_new, opt_pbpe%dm_total)

      opt_pbpe%norm_dm_1_diff = cp_fm_norm(opt_pbpe%dm_1_diff, 'F')
      opt_pbpe%norm_dm_2_diff = cp_fm_norm(opt_pbpe%dm_2_diff, 'F')
      opt_pbpe%norm_dm_total_diff = cp_fm_norm(opt_pbpe%dm_total_diff, 'F')

   END SUBROUTINE calculate_pbpe_dm_difference

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE check_pbpe_low_level_convergence(opt_pbpe)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_pbpe_low_level_convergence'

      IF (opt_pbpe%norm_dm_1_diff .LT. opt_pbpe%conv_rmse .AND. &
          opt_pbpe%norm_dm_2_diff .LT. opt_pbpe%conv_rmse .AND. &
          opt_pbpe%norm_dm_total_diff .LT. opt_pbpe%conv_rmse) THEN
         opt_pbpe%low_level_converged = .TRUE.
      ELSE
         opt_pbpe%low_level_converged = .FALSE.
      END IF

   END SUBROUTINE check_pbpe_low_level_convergence

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
! **************************************************************************************************
   SUBROUTINE check_pbpe_system(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_pbpe_system'

      INTEGER                                            :: handle, i_force_eval, n_force_eval, &
                                                            cluster_subsys_num, high_level_subsys_num, &
                                                            ref_subsys_num
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: naos, nmos
      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

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
      IF (cluster_subsys_num .NE. 2) THEN
         CPABORT("Cluster subsystem must be the second FORCE_EVAL.")
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

      CALL timestop(handle)
   
   END SUBROUTINE check_pbpe_system

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE copy_ks_to_pbpe(qs_env, opt_pbpe)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'copy_ks_to_pbpe'

      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ks_kp

      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks_kp)
      CALL copy_dbcsr_to_fm(matrix_ks_kp(1, 1)%matrix, opt_pbpe%ks_total)

      CALL cp_fm_to_fm_submat(opt_pbpe%ks_total, opt_pbpe%ks_1, opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(1), &
                              1, 1, 1, 1)
      CALL cp_fm_to_fm_submat(opt_pbpe%ks_total, opt_pbpe%ks_2, opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(2), &
                              opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, 1, 1)

      CALL cp_fm_to_fm_submat(opt_pbpe%ks_total, opt_pbpe%ks_12, opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(2), &
                              1, opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, 1, 1)
      CALL cp_fm_to_fm_submat(opt_pbpe%ks_total, opt_pbpe%ks_21, opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(1), &
                              opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, 1, 1, 1)

      END SUBROUTINE copy_ks_to_pbpe

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE create_and_init_pbpe_mo_sets(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_and_init_pbpe_mo_sets'

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

      ALLOCATE (opt_pbpe%mos_1(1))
      ALLOCATE (opt_pbpe%mos_2(1))

      NULLIFY (mo_coeff_1, mo_coeff_2)
      NULLIFY (fm_struct)

      ! Initialize mos_1
      CALL get_qs_env(qs_env=force_env%sub_force_env(1)%force_env%qs_env, mos=mos)
      CALL duplicate_mo_set(opt_pbpe%mos_1(1), mos(1))
      ! Copy mos_1(1)%mo_coeff to mo_coeff_1
      CALL get_mo_set(opt_pbpe%mos_1(1), mo_coeff=mo_coeff)
      ALLOCATE (mo_coeff_1)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(1), ncol_global=opt_pbpe%nmos(1))
      CALL cp_fm_create(mo_coeff_1, fm_struct)
      CALL cp_fm_set_all(mo_coeff_1, 0.0_dp)
      CALL cp_fm_to_fm_submat(mo_coeff, mo_coeff_1, opt_pbpe%naos(1), opt_pbpe%nmos(1), &
                              1, 1, 1, 1)
      ! Replace mos_1(1)%mo_coeff with mo_coeff_1
      CALL set_mo_set(opt_pbpe%mos_1(1), nao=opt_pbpe%naos_in_pbpe(1))
      CALL cp_fm_release(opt_pbpe%mos_1(1)%mo_coeff)
      DEALLOCATE (opt_pbpe%mos_1(1)%mo_coeff)
      opt_pbpe%mos_1(1)%mo_coeff => mo_coeff_1

      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      ! Initialize mos_2
      CALL get_qs_env(qs_env=force_env%sub_force_env(2)%force_env%qs_env, mos=mos)
      CALL duplicate_mo_set(opt_pbpe%mos_2(1), mos(1))
      ! Copy mos_2(1)%mo_coeff to mo_coeff_2
      CALL get_mo_set(opt_pbpe%mos_2(1), mo_coeff=mo_coeff)
      ALLOCATE (mo_coeff_2)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(2), ncol_global=opt_pbpe%nmos(2))
      CALL cp_fm_create(mo_coeff_2, fm_struct)
      CALL cp_fm_set_all(mo_coeff_2, 0.0_dp)
      CALL cp_fm_to_fm_submat(mo_coeff, mo_coeff_2, opt_pbpe%naos(2), opt_pbpe%nmos(2), &
                              1, 1, opt_pbpe%naos_in_pbpe(2) - opt_pbpe%naos(2) + 1, 1)
      ! Replace mos_2(1)%mo_coeff with mo_coeff_2
      CALL set_mo_set(opt_pbpe%mos_2(1), nao=opt_pbpe%naos_in_pbpe(2))
      CALL cp_fm_release(opt_pbpe%mos_2(1)%mo_coeff)
      DEALLOCATE (opt_pbpe%mos_2(1)%mo_coeff)
      opt_pbpe%mos_2(1)%mo_coeff => mo_coeff_2

      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_blacs_env_release(blacs_env)

      ! Associate mos pointers
      opt_pbpe%mos_cluster => opt_pbpe%mos_2
      opt_pbpe%mos_environment => opt_pbpe%mos_1

   END SUBROUTINE create_and_init_pbpe_mo_sets

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE create_pbpe_matrices(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_pbpe_matrices'

      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL get_qs_env(qs_env=force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env, &
                      para_env=para_env)

      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      NULLIFY (opt_pbpe%dm_1, opt_pbpe%dm_2, opt_pbpe%dm_total)
      NULLIFY (opt_pbpe%dm_1_new, opt_pbpe%dm_2_new, opt_pbpe%dm_total_new)
      NULLIFY (opt_pbpe%dm_1_beta, opt_pbpe%dm_2_beta, opt_pbpe%dm_total_beta)
      NULLIFY (opt_pbpe%dm_1_beta_new, opt_pbpe%dm_2_beta_new, opt_pbpe%dm_total_beta_new)
      NULLIFY (opt_pbpe%dm_1_diff, opt_pbpe%dm_2_diff, opt_pbpe%dm_total_diff)
      NULLIFY (opt_pbpe%dm_1_diff_beta, opt_pbpe%dm_2_diff_beta, opt_pbpe%dm_total_diff_beta)
      NULLIFY (opt_pbpe%dm_1_canonical, opt_pbpe%dm_2_canonical, opt_pbpe%dm_total_canonical)
      NULLIFY (opt_pbpe%dm_cluster, opt_pbpe%dm_environment)
      NULLIFY (opt_pbpe%dm_cluster_beta, opt_pbpe%dm_environment_beta)
      NULLIFY (opt_pbpe%dm_cluster_diff, opt_pbpe%dm_environment_diff)
      NULLIFY (opt_pbpe%dm_cluster_diff_beta, opt_pbpe%dm_environment_diff_beta)
      NULLIFY (opt_pbpe%dm_cluster_canonical, opt_pbpe%dm_environment_canonical)
      NULLIFY (opt_pbpe%s_total, opt_pbpe%ks_total)
      NULLIFY (opt_pbpe%s_1, opt_pbpe%s_2)
      NULLIFY (opt_pbpe%ks_1, opt_pbpe%ks_2)
      NULLIFY (opt_pbpe%ks_1_in_2, opt_pbpe%ks_2_in_1)
      NULLIFY (opt_pbpe%s_12, opt_pbpe%s_21)
      NULLIFY (opt_pbpe%ks_12, opt_pbpe%ks_21)
      NULLIFY (opt_pbpe%work1_1, opt_pbpe%work1_2)
      NULLIFY (opt_pbpe%work2_1, opt_pbpe%work2_2)
      NULLIFY (opt_pbpe%ortho_1, opt_pbpe%ortho_2)
      NULLIFY (opt_pbpe%projection_1, opt_pbpe%projection_2)
      NULLIFY (opt_pbpe%projection_environment) 

      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(3), ncol_global=opt_pbpe%naos_in_pbpe(3))
      ALLOCATE (opt_pbpe%dm_total, opt_pbpe%dm_total_new, opt_pbpe%dm_total_diff, opt_pbpe%dm_total_canonical, &
                opt_pbpe%s_total, opt_pbpe%ks_total)
      CALL cp_fm_create(opt_pbpe%dm_total, fm_struct, name="dm_total")
      CALL cp_fm_create(opt_pbpe%dm_total_new, fm_struct, name="dm_total_new")
      CALL cp_fm_create(opt_pbpe%dm_total_diff, fm_struct, name="dm_total_diff")
      CALL cp_fm_create(opt_pbpe%dm_total_canonical, fm_struct, name="dm_total_canonical")
      CALL cp_fm_create(opt_pbpe%s_total, fm_struct, name="s_total")
      CALL cp_fm_create(opt_pbpe%ks_total, fm_struct, name="KOHN-SHAM MATRIX")

      CALL cp_fm_set_all(opt_pbpe%dm_total, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_total_new, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_total_diff, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_total_canonical, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%s_total, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_total, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_pbpe%open_shell_embed) THEN
         ALLOCATE (opt_pbpe%dm_total_beta, opt_pbpe%dm_total_beta_new, opt_pbpe%dm_total_diff_beta)
         CALL cp_fm_create(opt_pbpe%dm_total_beta, fm_struct, name="dm_total_beta")
         CALL cp_fm_create(opt_pbpe%dm_total_beta_new, fm_struct, name="dm_total_beta_new")
         CALL cp_fm_create(opt_pbpe%dm_total_diff_beta, fm_struct, name="dm_total_diff_beta")

         CALL cp_fm_set_all(opt_pbpe%dm_total_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_total_beta_new, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_total_diff_beta, 0.0_dp)
      END IF

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(1), ncol_global=opt_pbpe%naos_in_pbpe(1))
      ALLOCATE (opt_pbpe%dm_1, opt_pbpe%dm_1_new, opt_pbpe%dm_1_diff, opt_pbpe%dm_1_canonical, &
                opt_pbpe%s_1, opt_pbpe%ks_1, opt_pbpe%ks_1_in_2, opt_pbpe%work1_1, opt_pbpe%work2_1, &
                opt_pbpe%ortho_1, opt_pbpe%projection_2, opt_pbpe%fm_tmp3_1)
      CALL cp_fm_create(opt_pbpe%dm_1, fm_struct, name="DENSITY MATRIX OF SUBSYSTEM 1")
      CALL cp_fm_create(opt_pbpe%dm_1_new, fm_struct, name="dm_1_new")
      CALL cp_fm_create(opt_pbpe%dm_1_diff, fm_struct, name="dm_1_diff")
      CALL cp_fm_create(opt_pbpe%dm_1_canonical, fm_struct, name="dm_1_canonical")
      CALL cp_fm_create(opt_pbpe%s_1, fm_struct, name="OVERLAP MATRIX OF SUBSYSTEM 1")
      CALL cp_fm_create(opt_pbpe%ks_1, fm_struct, name="KOHN-SHAM MATRIX OF SUBSYSTEM 1")
      CALL cp_fm_create(opt_pbpe%ks_1_in_2, fm_struct, name="EMBEDDED KOHN-SHAM MATRIX OF SUBSYSTEM 1")
      CALL cp_fm_create(opt_pbpe%work1_1, fm_struct, name="work1_1")
      CALL cp_fm_create(opt_pbpe%work2_1, fm_struct, name="work2_1")
      CALL cp_fm_create(opt_pbpe%ortho_1, fm_struct, name="ORTHO OF SUBSYSTEM 1")
      CALL cp_fm_create(opt_pbpe%projection_2, fm_struct, name="projection_2")
      CALL cp_fm_create(opt_pbpe%fm_tmp3_1, fm_struct, name="(ks12-e*s12)*d2*s21")

      CALL cp_fm_set_all(opt_pbpe%dm_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_1_new, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_1_diff, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%s_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_1_in_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work1_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work2_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ortho_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%projection_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%fm_tmp3_1, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_pbpe%open_shell_embed) THEN
         ALLOCATE (opt_pbpe%dm_1_beta, opt_pbpe%dm_1_beta_new, opt_pbpe%dm_1_diff_beta)
         CALL cp_fm_create(opt_pbpe%dm_1_beta, fm_struct, name="dm_1_beta")
         CALL cp_fm_create(opt_pbpe%dm_1_beta_new, fm_struct, name="dm_1_beta_new")
         CALL cp_fm_create(opt_pbpe%dm_1_diff_beta, fm_struct, name="dm_1_diff_beta")

         CALL cp_fm_set_all(opt_pbpe%dm_1_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_1_beta_new, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_1_diff_beta, 0.0_dp)
      END IF

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(2), ncol_global=opt_pbpe%naos_in_pbpe(2))
      ALLOCATE (opt_pbpe%dm_2, opt_pbpe%dm_2_new, opt_pbpe%dm_2_diff, opt_pbpe%dm_2_canonical, &
                opt_pbpe%s_2, opt_pbpe%ks_2, opt_pbpe%ks_2_in_1, opt_pbpe%work1_2, opt_pbpe%work2_2, &
                opt_pbpe%ortho_2, opt_pbpe%projection_1, opt_pbpe%fm_tmp3_2)
      CALL cp_fm_create(opt_pbpe%dm_2, fm_struct, name="DENSITY MATRIX OF SUBSYSTEM 2")
      CALL cp_fm_create(opt_pbpe%dm_2_new, fm_struct, name="dm_2_new")
      CALL cp_fm_create(opt_pbpe%dm_2_diff, fm_struct, name="dm_2_diff")
      CALL cp_fm_create(opt_pbpe%dm_2_canonical, fm_struct, name="dm_2_canonical")
      CALL cp_fm_create(opt_pbpe%s_2, fm_struct, name="OVERLAP MATRIX OF SUBSYSTEM 2")
      CALL cp_fm_create(opt_pbpe%ks_2, fm_struct, name="KOHN-SHAM MATRIX OF SUBSYSTEM 2")
      CALL cp_fm_create(opt_pbpe%ks_2_in_1, fm_struct, name="EMBEDDED KOHN-SHAM MATRIX OF SUBSYSTEM 2")
      CALL cp_fm_create(opt_pbpe%work1_2, fm_struct, name="work1_2")
      CALL cp_fm_create(opt_pbpe%work2_2, fm_struct, name="work2_2")
      CALL cp_fm_create(opt_pbpe%ortho_2, fm_struct, name="ORTHO OF SUBSYSTEM 2")
      CALL cp_fm_create(opt_pbpe%projection_1, fm_struct, name="projection_1")
      CALL cp_fm_create(opt_pbpe%fm_tmp3_2, fm_struct, name="(ks21-e*s21)*d1*s12")

      CALL cp_fm_set_all(opt_pbpe%dm_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_2_new, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%dm_2_diff, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%s_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_2_in_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work1_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%work2_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ortho_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%projection_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%fm_tmp3_2, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_pbpe%open_shell_embed) THEN
         ALLOCATE (opt_pbpe%dm_2_beta, opt_pbpe%dm_2_beta_new, opt_pbpe%dm_2_diff_beta)
         CALL cp_fm_create(opt_pbpe%dm_2_beta, fm_struct, name="dm_2_beta")
         CALL cp_fm_create(opt_pbpe%dm_2_beta_new, fm_struct, name="dm_2_beta")
         CALL cp_fm_create(opt_pbpe%dm_2_diff_beta, fm_struct, name="dm_2_diff_beta")

         CALL cp_fm_set_all(opt_pbpe%dm_2_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_2_beta_new, 0.0_dp)
         CALL cp_fm_set_all(opt_pbpe%dm_2_diff_beta, 0.0_dp)
      END IF

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(1), ncol_global=opt_pbpe%naos_in_pbpe(2))
      ALLOCATE (opt_pbpe%s_12, opt_pbpe%ks_12, opt_pbpe%fm_tmp1_1, opt_pbpe%fm_tmp2_1)
      CALL cp_fm_create(opt_pbpe%s_12, fm_struct, name="s_12")
      CALL cp_fm_create(opt_pbpe%ks_12, fm_struct, name="ks_12")
      CALL cp_fm_create(opt_pbpe%fm_tmp1_1, fm_struct, name="ks12-e*s12")
      CALL cp_fm_create(opt_pbpe%fm_tmp2_1, fm_struct, name="(ks12-e*s12)*d2")

      CALL cp_fm_set_all(opt_pbpe%s_12, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_12, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%fm_tmp1_1, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%fm_tmp2_1, 0.0_dp)

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=opt_pbpe%naos_in_pbpe(2), ncol_global=opt_pbpe%naos_in_pbpe(1))
      ALLOCATE (opt_pbpe%s_21, opt_pbpe%ks_21, opt_pbpe%fm_tmp1_2, opt_pbpe%fm_tmp2_2)
      CALL cp_fm_create(opt_pbpe%s_21, fm_struct, name="s_21")
      CALL cp_fm_create(opt_pbpe%ks_21, fm_struct, name="ks_21")
      CALL cp_fm_create(opt_pbpe%fm_tmp1_2, fm_struct, name="ks2-e*s21")
      CALL cp_fm_create(opt_pbpe%fm_tmp2_2, fm_struct, name="(ks2-e*s21)*d1")

      CALL cp_fm_set_all(opt_pbpe%s_21, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%ks_21, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%fm_tmp1_2, 0.0_dp)
      CALL cp_fm_set_all(opt_pbpe%fm_tmp2_2, 0.0_dp)

      ! Release structure
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)

      ! Release context
      CALL cp_blacs_env_release(blacs_env)

      ! Associate matrix pointers
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

   END SUBROUTINE create_pbpe_matrices

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE init_pbpe_matrices(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_pbpe_matrices'

      INTEGER                                            :: n_dep_1, n_dep_2
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_array_2d
      TYPE(qs_rho_type), POINTER                         :: rho

      ! Initialize density matrices
      CALL build_pbpe_subsys_dm(opt_pbpe, subsys_num=1, density_matrix=opt_pbpe%dm_1)
      CALL build_pbpe_subsys_dm(opt_pbpe, subsys_num=2, density_matrix=opt_pbpe%dm_2)
      CALL build_pbpe_ref_subsys_dm(opt_pbpe, density_matrix=opt_pbpe%dm_total, use_new_subsys_dm=.FALSE.)

      ! Initialize overlap matrices
      CALL get_qs_env(force_env%sub_force_env(3)%force_env%qs_env, matrix_s_kp=matrix_array_2d)
      CALL copy_dbcsr_to_fm(matrix_array_2d(1, 1)%matrix, opt_pbpe%s_total)
      CALL cp_fm_to_fm_submat(opt_pbpe%s_total, opt_pbpe%s_1, opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(1), &
                              1, 1, 1, 1)
      CALL cp_fm_to_fm_submat(opt_pbpe%s_total, opt_pbpe%s_2, opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(2), &
                              opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, 1, 1)
      CALL cp_fm_to_fm_submat(opt_pbpe%s_total, opt_pbpe%s_12, opt_pbpe%naos_in_pbpe(1), opt_pbpe%naos_in_pbpe(2), &
                              1, opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, 1, 1)
      CALL cp_fm_to_fm_submat(opt_pbpe%s_total, opt_pbpe%s_21, opt_pbpe%naos_in_pbpe(2), opt_pbpe%naos_in_pbpe(1), &
                              opt_pbpe%naos_in_pbpe(3) - opt_pbpe%naos_in_pbpe(2) + 1, 1, 1, 1)

      ! Prepare for diagonalization
      n_dep_1 = 0
      n_dep_2 = 0
      CALL cp_fm_to_fm(opt_pbpe%s_1, opt_pbpe%ortho_1)
      CALL cp_fm_power(opt_pbpe%ortho_1, opt_pbpe%work2_1, -0.5_dp, 1.0E-5_dp, n_dep_1)
      CALL cp_fm_to_fm(opt_pbpe%s_2, opt_pbpe%ortho_2)
      CALL cp_fm_power(opt_pbpe%ortho_2, opt_pbpe%work2_2, -0.5_dp, 1.0E-5_dp, n_dep_2)

   END SUBROUTINE init_pbpe_matrices

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE pbpe_do_diagonalization(opt_pbpe, subsys_num)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      INTEGER, INTENT(IN)                                :: subsys_num

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pbpe_do_diagonalization'

      CPASSERT((subsys_num .EQ. 1) .OR. (subsys_num .EQ. 2))
      IF (subsys_num .EQ. 1) THEN
         CALL cp_fm_to_fm(opt_pbpe%ks_1_in_2, opt_pbpe%work1_1)
         CALL eigensolver_symm(matrix_ks_fm=opt_pbpe%work1_1, mo_set=opt_pbpe%mos_1(1), &
                               ortho=opt_pbpe%ortho_1, work=opt_pbpe%work2_1, &
                               do_level_shift=.FALSE., level_shift=0.0_dp, &
                               use_jacobi=.FALSE., jacobi_threshold=0.0_dp)
      ELSE IF (subsys_num .EQ. 2) THEN
         CALL cp_fm_to_fm(opt_pbpe%ks_2_in_1, opt_pbpe%work1_2)
         CALL eigensolver_symm(matrix_ks_fm=opt_pbpe%work1_2, mo_set=opt_pbpe%mos_2(1), &
                               ortho=opt_pbpe%ortho_2, work=opt_pbpe%work2_2, &
                               do_level_shift=.FALSE., level_shift=0.0_dp, &
                               use_jacobi=.FALSE., jacobi_threshold=0.0_dp)
      END IF

   END SUBROUTINE pbpe_do_diagonalization

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE pbpe_initialization_print(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pbpe_initial_print'

      INTEGER                                            :: after, iw
      LOGICAL                                            :: omit_headers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), POINTER                 :: qs_env_1, qs_env_2, qs_env_ref
      TYPE(section_vals_type), POINTER                   :: input

      logger => cp_get_default_logger()

      qs_env_ref => force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env
      CALL get_qs_env(qs_env_ref, input=input, para_env=para_env)

      IF (opt_pbpe%abs_localization) THEN
         qs_env_1 => force_env%sub_force_env(1)%force_env%qs_env
         qs_env_2 => force_env%sub_force_env(2)%force_env%qs_env
      ELSE
         qs_env_1 => force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env
         qs_env_2 => force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env
      END IF

      CALL section_vals_val_get(input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES%NDIGITS", i_val=after)
      CALL section_vals_val_get(input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      after = MIN(MAX(after, 1), 16)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/ORTHO"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/ORTHO", extension=".Log")
         CALL write_fm_with_basis_info(opt_pbpe%ortho_1, 4, after, qs_env=qs_env_1, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL write_fm_with_basis_info(opt_pbpe%ortho_2, 4, after, qs_env=qs_env_2, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL cp_print_key_finished_output(iw, logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/ORTHO")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/OVERLAP"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/OVERLAP", extension=".Log")
         CALL write_fm_with_basis_info(opt_pbpe%s_1, 4, after, qs_env=qs_env_1, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL write_fm_with_basis_info(opt_pbpe%s_2, 4, after, qs_env=qs_env_2, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL cp_print_key_finished_output(iw, logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/OVERLAP")
      END IF

   END SUBROUTINE pbpe_initialization_print

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE pbpe_low_level_loop_print(force_env, opt_pbpe)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pbpe_low_level_loop_print'

      INTEGER                                            :: after, iw
      LOGICAL                                            :: omit_headers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), POINTER                 :: qs_env_1, qs_env_2, qs_env_ref
      TYPE(section_vals_type), POINTER                   :: input

      logger => cp_get_default_logger()

      qs_env_ref => force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env
      CALL get_qs_env(qs_env_ref, input=input, para_env=para_env)

      IF (opt_pbpe%abs_localization) THEN
         qs_env_1 => force_env%sub_force_env(1)%force_env%qs_env
         qs_env_2 => force_env%sub_force_env(2)%force_env%qs_env
      ELSE
         qs_env_1 => force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env
         qs_env_2 => force_env%sub_force_env(opt_pbpe%ref_subsys_num)%force_env%qs_env
      END IF

      CALL section_vals_val_get(input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES%NDIGITS", i_val=after)
      CALL section_vals_val_get(input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      after = MIN(MAX(after, 1), 16)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/DENSITY"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/DENSITY", extension=".Log")
         CALL write_fm_with_basis_info(opt_pbpe%dm_1, 4, after, qs_env=qs_env_1, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL write_fm_with_basis_info(opt_pbpe%dm_2, 4, after, qs_env=qs_env_2, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL cp_print_key_finished_output(iw, logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/DENSITY")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/EMBEDDED_KOHN_SHAM_MATRICES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/EMBEDDED_KOHN_SHAM_MATRICES", extension=".Log")
         CALL write_fm_with_basis_info(opt_pbpe%ks_1_in_2, 4, after, qs_env=qs_env_1, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL write_fm_with_basis_info(opt_pbpe%ks_2_in_1, 4, after, qs_env=qs_env_2, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL cp_print_key_finished_output(iw, logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/EMBEDDED_KOHN_SHAM_MATRICES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", extension=".Log")
         CALL write_fm_with_basis_info(opt_pbpe%ks_total, 4, after, qs_env=qs_env_ref, &
                                       para_env=para_env, output_unit=iw, omit_headers=omit_headers)
         CALL cp_print_key_finished_output(iw, logger, input, "DFT%QS%OPT_PBPE%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX")
      END IF

   END SUBROUTINE pbpe_low_level_loop_print

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE pbpe_summary_print(opt_pbpe, output_unit)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      INTEGER                                            :: output_unit
   
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pbpe_summary_print'


   END SUBROUTINE pbpe_summary_print

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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prepare_pbpe_opt'

      INTEGER                                            :: handle, i_force_eval, n_force_eval
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

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

      opt_pbpe%environment_subsys_num = 1

      DO i_force_eval = 1, n_force_eval
         opt_pbpe%naos_in_pbpe(i_force_eval) = opt_pbpe%naos(i_force_eval) 
      END DO

      ! modify naos_in_pbpe if absolute localization is not used
      IF (.NOT. opt_pbpe%abs_localization) THEN
         opt_pbpe%naos_in_pbpe(opt_pbpe%cluster_subsys_num) = opt_pbpe%naos_in_pbpe(opt_pbpe%ref_subsys_num)
         opt_pbpe%naos_in_pbpe(opt_pbpe%environment_subsys_num) = opt_pbpe%naos_in_pbpe(opt_pbpe%ref_subsys_num)
         opt_pbpe%naos_in_pbpe(opt_pbpe%high_level_subsys_num) = opt_pbpe%naos_in_pbpe(opt_pbpe%ref_subsys_num)
      END IF

      ! Create and initialize MO sets
      CALL create_and_init_pbpe_mo_sets(force_env, opt_pbpe)

      ! Create and initialize cp_fm matrices
      CALL create_pbpe_matrices(force_env, opt_pbpe)
      CALL init_pbpe_matrices(force_env, opt_pbpe)

      ! Get Fermi energy
      CALL get_mo_set(opt_pbpe%mos_1(1), mu=opt_pbpe%homo_energy_1)
      CALL get_mo_set(opt_pbpe%mos_2(1), mu=opt_pbpe%homo_energy_2)
      IF (opt_pbpe%homo_energy_1 .GT. opt_pbpe%homo_energy_2) THEN
         opt_pbpe%fermi_energy = opt_pbpe%homo_energy_1
      ELSE
         opt_pbpe%fermi_energy = opt_pbpe%homo_energy_2
      END IF
      IF (opt_pbpe%fermi_energy .LT. 0) THEN
         opt_pbpe%fermi_energy = 0
      END IF

      CALL timestop(handle)

   END SUBROUTINE prepare_pbpe_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE read_opt_pbpe_section(opt_pbpe, opt_pbpe_section)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe
      TYPE(section_vals_type), POINTER                   :: opt_pbpe_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_opt_pbpe_section'

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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'release_pbpe_matrices'

      NULLIFY (opt_pbpe%dm_cluster)
      NULLIFY (opt_pbpe%dm_environment)
      NULLIFY (opt_pbpe%dm_cluster_diff)
      NULLIFY (opt_pbpe%dm_environment_diff)
      NULLIFY (opt_pbpe%dm_cluster_canonical)
      NULLIFY (opt_pbpe%dm_environment_canonical)

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
      IF (ASSOCIATED(opt_pbpe%s_total)) THEN
         CALL cp_fm_release(opt_pbpe%s_total)
         DEALLOCATE (opt_pbpe%s_total)
         NULLIFY (opt_pbpe%s_total)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_total)) THEN
         CALL cp_fm_release(opt_pbpe%ks_total)
         DEALLOCATE (opt_pbpe%ks_total)
         NULLIFY (opt_pbpe%ks_total)
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
      IF (ASSOCIATED(opt_pbpe%dm_1_canonical)) THEN
         CALL cp_fm_release(opt_pbpe%dm_1_canonical)
         DEALLOCATE (opt_pbpe%dm_1_canonical)
         NULLIFY (opt_pbpe%dm_1_canonical)
      END IF
      IF (ASSOCIATED(opt_pbpe%s_1)) THEN
         CALL cp_fm_release(opt_pbpe%s_1)
         DEALLOCATE (opt_pbpe%s_1)
         NULLIFY (opt_pbpe%s_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_1)) THEN
         CALL cp_fm_release(opt_pbpe%ks_1)
         DEALLOCATE (opt_pbpe%ks_1)
         NULLIFY (opt_pbpe%ks_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_1_in_2)) THEN
         CALL cp_fm_release(opt_pbpe%ks_1_in_2)
         DEALLOCATE (opt_pbpe%ks_1_in_2)
         NULLIFY (opt_pbpe%ks_1_in_2)
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
      IF (ASSOCIATED(opt_pbpe%projection_2)) THEN
         CALL cp_fm_release(opt_pbpe%projection_2)
         DEALLOCATE (opt_pbpe%projection_2)
         NULLIFY (opt_pbpe%projection_2)
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
      IF (ASSOCIATED(opt_pbpe%dm_2_canonical)) THEN
         CALL cp_fm_release(opt_pbpe%dm_2_canonical)
         DEALLOCATE (opt_pbpe%dm_2_canonical)
         NULLIFY (opt_pbpe%dm_2_canonical)
      END IF
      IF (ASSOCIATED(opt_pbpe%s_2)) THEN
         CALL cp_fm_release(opt_pbpe%s_2)
         DEALLOCATE (opt_pbpe%s_2)
         NULLIFY (opt_pbpe%s_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_2)) THEN
         CALL cp_fm_release(opt_pbpe%ks_2)
         DEALLOCATE (opt_pbpe%ks_2)
         NULLIFY (opt_pbpe%ks_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_2_in_1)) THEN
         CALL cp_fm_release(opt_pbpe%ks_2_in_1)
         DEALLOCATE (opt_pbpe%ks_2_in_1)
         NULLIFY (opt_pbpe%ks_2_in_1)
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
      IF (ASSOCIATED(opt_pbpe%projection_1)) THEN
         CALL cp_fm_release(opt_pbpe%projection_1)
         DEALLOCATE (opt_pbpe%projection_1)
         NULLIFY (opt_pbpe%projection_1)
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

      IF (ASSOCIATED(opt_pbpe%s_12)) THEN
         CALL cp_fm_release(opt_pbpe%s_12)
         DEALLOCATE (opt_pbpe%s_12)
         NULLIFY (opt_pbpe%s_12)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_12)) THEN
         CALL cp_fm_release(opt_pbpe%ks_12)
         DEALLOCATE (opt_pbpe%ks_12)
         NULLIFY (opt_pbpe%ks_12)
      END IF

      IF (ASSOCIATED(opt_pbpe%s_21)) THEN
         CALL cp_fm_release(opt_pbpe%s_21)
         DEALLOCATE (opt_pbpe%s_21)
         NULLIFY (opt_pbpe%s_21)
      END IF
      IF (ASSOCIATED(opt_pbpe%ks_21)) THEN
         CALL cp_fm_release(opt_pbpe%ks_21)
         DEALLOCATE (opt_pbpe%ks_21)
         NULLIFY (opt_pbpe%ks_21)
      END IF

      IF (ASSOCIATED(opt_pbpe%fm_tmp1_1)) THEN
         CALL cp_fm_release(opt_pbpe%fm_tmp1_1)
         DEALLOCATE (opt_pbpe%fm_tmp1_1)
         NULLIFY (opt_pbpe%fm_tmp1_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%fm_tmp2_1)) THEN
         CALL cp_fm_release(opt_pbpe%fm_tmp2_1)
         DEALLOCATE (opt_pbpe%fm_tmp2_1)
         NULLIFY (opt_pbpe%fm_tmp2_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%fm_tmp3_1)) THEN
         CALL cp_fm_release(opt_pbpe%fm_tmp3_1)
         DEALLOCATE (opt_pbpe%fm_tmp3_1)
         NULLIFY (opt_pbpe%fm_tmp3_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%fm_tmp1_2)) THEN
         CALL cp_fm_release(opt_pbpe%fm_tmp1_2)
         DEALLOCATE (opt_pbpe%fm_tmp1_2)
         NULLIFY (opt_pbpe%fm_tmp1_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%fm_tmp2_2)) THEN
         CALL cp_fm_release(opt_pbpe%fm_tmp2_2)
         DEALLOCATE (opt_pbpe%fm_tmp2_2)
         NULLIFY (opt_pbpe%fm_tmp2_2)
      END IF
      IF (ASSOCIATED(opt_pbpe%fm_tmp3_2)) THEN
         CALL cp_fm_release(opt_pbpe%fm_tmp3_2)
         DEALLOCATE (opt_pbpe%fm_tmp3_2)
         NULLIFY (opt_pbpe%fm_tmp3_2)
      END IF

   END SUBROUTINE release_pbpe_matrices

! **************************************************************************************************
!> \brief ...
!> \param opt_pbpe ...
! **************************************************************************************************
   SUBROUTINE release_pbpe_mo_sets(opt_pbpe)
      TYPE(opt_pbpe_pot_type)                            :: opt_pbpe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'release_pbpe_mo_sets'

      NULLIFY (opt_pbpe%mos_cluster)
      NULLIFY (opt_pbpe%mos_environment)

      IF (ASSOCIATED(opt_pbpe%mos_1)) THEN
         CALL deallocate_mo_set(opt_pbpe%mos_1(1))
         DEALLOCATE (opt_pbpe%mos_1)
         NULLIFY (opt_pbpe%mos_1)
      END IF
      IF (ASSOCIATED(opt_pbpe%mos_2)) THEN
         CALL deallocate_mo_set(opt_pbpe%mos_2(1))
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'release_pbpe_opt'

      CALL release_pbpe_matrices(opt_pbpe)

      CALL release_pbpe_mo_sets(opt_pbpe)

      DEALLOCATE (opt_pbpe%naos)
      DEALLOCATE (opt_pbpe%naos_in_pbpe)
      DEALLOCATE (opt_pbpe%nmos)

   END SUBROUTINE release_pbpe_opt

! **************************************************************************************************
!> \brief ...
!> \param density_matrix ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE update_pbpe_density_and_ks_qs_env(density_matrix, qs_env)
   TYPE(cp_fm_type), INTENT(IN)                          :: density_matrix
   TYPE(qs_environment_type), POINTER                    :: qs_env

   CHARACTER(LEN=*), PARAMETER                           :: routineN = 'update_pbpe_density_and_ks_qs_env'

   TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER           :: rho_ao_kp
   TYPE(qs_ks_env_type), POINTER                         :: ks_env
   TYPE(qs_rho_type), POINTER                            :: rho
   TYPE(qs_scf_env_type), POINTER                        :: scf_env 

   CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho, scf_env=scf_env)
   CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
   CALL copy_fm_to_dbcsr(density_matrix, rho_ao_kp(1, 1)%matrix)
   CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, .FALSE.)
   CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., calculate_forces=.FALSE., print_active=.TRUE.)

   END SUBROUTINE update_pbpe_density_and_ks_qs_env

END MODULE optimize_pbpe_potential