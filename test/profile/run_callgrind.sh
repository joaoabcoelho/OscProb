#!/bin/bash

if [ $# -eq 0 ]; then
  MODEL=Fast
else
  MODEL=$1
fi

make

echo "Running profiler over PMNS_$MODEL"

valgrind --tool=callgrind --callgrind-out-file=callgrind.out -q ./bin/stress $MODEL

echo "Open profiler output"

if command -v kcachegrind &> /dev/null
then
  kcachegrind callgrind.out
elif command -v callgrind_annotate &> /dev/null
then
  callgrind_annotate callgrind.out
else
  echo "Neither kcachegrind nor callgrind_annotate executables found"
  echo "Cannot open callgrind output"
fi
