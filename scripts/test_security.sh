#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
test_directory=$(mktemp -d)
trap 'rm -rf "$test_directory"' EXIT

flags=(-std=c++14 -O1 -g -DNDEBUG -fno-omit-frame-pointer
       -fsanitize=address,undefined -fno-sanitize-recover=all -Isrc)
"${CXX:-g++}" "${flags[@]}" tests/security_inputs.cpp -o "$test_directory/security_inputs"
"$test_directory/security_inputs" "$test_directory"

expect_rejection() {
    local status=0
    "$@" > "$test_directory/rejected.log" 2>&1 || status=$?
    if [[ "$status" != 1 ]] || ! grep -q '^Error:' "$test_directory/rejected.log"; then
        cat "$test_directory/rejected.log"
        printf 'Expected a controlled error (exit 1), got %s\n' "$status" >&2
        exit 1
    fi
}

printf 'ATOM 1 CA GLY A 1 2\n' > "$test_directory/malformed.pdb"
for program in execute_sh_ilo execute_shibuya_method; do
    "${CXX:-g++}" "${flags[@]}" "src/$program.cpp" -o "$test_directory/$program"
    binary="$test_directory/$program"
    pair=(src/pdb3hvp.pdb src/pdb4hvp.pdb A A)
    "$binary" "${pair[@]}" 2
    for hinges in -1 0 999999999999999999999 1000000 invalid 2junk ''; do
        expect_rejection "$binary" "${pair[@]}" "$hinges"
    done
    expect_rejection "$binary" src/pdb3hvp.pdb src/pdb4hvp.pdb '' A 2
    expect_rejection "$binary" src/pdb3hvp.pdb src/pdb4hvp.pdb missing A 2
    expect_rejection "$binary" "$test_directory/missing.pdb" src/pdb4hvp.pdb A A 2
    expect_rejection "$binary" "$test_directory/malformed.pdb" src/pdb4hvp.pdb A A 2
done
printf 'CLI security regression tests passed\n'
