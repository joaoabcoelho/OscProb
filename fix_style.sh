#!/bin/bash

SEARCHARG=(src inc \
          ! -name "LinkDef.h" \
          \( -name "*.cxx" -o -name "*.h" \))

if [ "$1" == "test" ]; then
  find ${SEARCHARG[@]} \
       -exec clang-format -n --Werror {} +
else
  find ${SEARCHARG[@]} \
       -exec clang-format -i --verbose {} +
fi
