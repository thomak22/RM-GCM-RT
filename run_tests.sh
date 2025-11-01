#!/bin/bash
cp params.i_hires params.i
make clean
make

for dir in tests*/; do
    if [[ "$dir" != *hires/ ]]; then
        continue
    fi
    if [[ -d "$dir" ]]; then
        files=()
        while IFS= read -r -d '' f; do
            files+=("$(basename "$f")")
        done < <(find "$dir" -maxdepth 1 -type f -print0)
    else
        files=()
    fi
    echo "Found test files in $dir: ${files[*]}"
    # cp "$dir/params.i" .
    for testfile in "${files[@]}"; do
        if [[ "$testfile" == "params.i" ]]; then
            continue
        fi
        echo "Running test with input file: $testfile in $dir"
        rm -r Planet_Run
        mkdir Planet_Run
        cp 1DRT Planet_Run
        cp "$dir/$testfile" Planet_Run/fort.7
        cd Planet_Run
        ./1DRT
        cd ..
        rm -r "Planet_Run_$testfile"
        mv Planet_Run "Planet_Run_$testfile"
    done
done

cp params.i_lores params.i
make clean
make

for dir in tests*/; do
    if [[ "$dir" == *hires/ ]]; then
        continue
    fi
    if [[ -d "$dir" ]]; then
        files=()
        while IFS= read -r -d '' f; do
            files+=("$(basename "$f")")
        done < <(find "$dir" -maxdepth 1 -type f -print0)
    else
        files=()
    fi
    echo "Found test files in $dir: ${files[*]}"
    # cp "$dir/params.i" .
    for testfile in "${files[@]}"; do
        if [[ "$testfile" == "params.i" ]]; then
            continue
        fi
        echo "Running test with input file: $testfile in $dir"
        rm -r Planet_Run
        mkdir Planet_Run
        cp 1DRT Planet_Run
        cp "$dir/$testfile" Planet_Run/fort.7
        cd Planet_Run
        ./1DRT
        cd ..
        rm -r "Planet_Run_$testfile"
        mv Planet_Run "Planet_Run_$testfile"
    done
done

# dir="testshires"
# if [[ -d "$dir" ]]; then
#     files=()
#     while IFS= read -r -d '' f; do
#         files+=("$(basename "$f")")
#     done < <(find "$dir" -maxdepth 1 -type f -print0)
# else
#     files=()
# fi
# echo "Found test files: ${files[*]}"
# cp "$dir/params.i" .
# make clean
# make
# for testfile in "${files[@]}"; do
#     if [[ "$testfile" == "params.i" ]]; then
#         continue
#     fi
#     echo "Running test with input file: $testfile"
#     rm -r Planet_Run
#     mkdir Planet_Run
#     cp 1DRT Planet_Run
#     cp "$dir/$testfile" Planet_Run/fort.7
#     cd Planet_Run
#     ./1DRT
#     cd ..
#     rm -r "Planet_Run_$testfile"
#     mv Planet_Run "Planet_Run_$testfile"
# done