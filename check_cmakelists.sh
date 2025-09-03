#!/usr/bin/env bash

src_dir=$(dirname "$(readlink -f "$0")")

function check() {
    dir="$src_dir/$1"
    extension="$2"
    files=()
    while IFS= read -r -d $'\0'; do
        files+=("$REPLY")
    done < <(find "$dir" -name "*$extension" -print0)

    in_cmakelists=$(find "$src_dir" -name "CMakeLists.txt" -exec grep "$extension" {} +)

    is_bad=0
    for filename in "${files[@]}"; do
        base_filename=$(basename "$filename")
        [[ $in_cmakelists == *"${base_filename}"* ]] || {
            echo "File $filename not found in any CMakeLists.txt"
            is_bad=1
        }
    done

    return "$is_bad"

}

check "." ".cxx" || exit 1
check "inc" ".h" || exit 1
