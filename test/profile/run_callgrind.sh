#!/bin/bash

if [ $# -eq 0 ]; then
  MODEL=Fast
else
  MODEL=$1
fi

echo "Running profiler over PMNS_$MODEL"

valgrind --tool=callgrind --callgrind-out-file=callgrind.out -q ./bin/stress $MODEL

echo "Open profiler output"

kcachegrind callgrind.out
