#!/bin/bash -e

# author: Ole Schuett

echo -e "\n========== Generating Doxygen =========="
cd /opt/cp2k
if ! make doxygen &> make.out; then
  echo -e "failed.\n\n"
  tail -n 100 make.out
  mkdir -p /workspace/artifacts/
  cp make.out /workspace/artifacts/
  echo -e "\nSummary: Doxygen generation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

mkdir -p /workspace/artifacts
mv ./doxygen/html /workspace/artifacts/doxygen

echo "Summary: Doxygen generation works fine."
echo "Status: OK"

#EOF
