#!/usr/bin/env bash
# локальное окружение Valgrind для CachyOS
set -euo pipefail

if [[ $(uname -s) != Linux || $(uname -m) != x86_64 ]]; then
    echo 'This runtime requires Linux x86_64.' >&2
    exit 1
fi

build_dir=${1:?Usage: setup_valgrind.sh BUILD_DIR}
runtime="$build_dir/valgrind-glibc"
mkdir -p "$runtime"

for package in glibc valgrind libstdc++ libgcc; do
    repository=core
    archive="$build_dir/$package.pkg.tar.zst"
    if [[ $package == valgrind ]]; then
        repository=extra
    elif [[ $package == glibc ]]; then
        archive="$build_dir/glibc-valgrind.pkg.tar.zst"
    fi
    if [[ ! -s $archive ]]; then
        curl --fail --location --retry 2 --max-time 120 \
            "https://archlinux.org/packages/$repository/x86_64/$package/download/" \
            --output "$archive.part"
        mv -- "$archive.part" "$archive"
    fi
    bsdtar -xf "$archive" -C "$runtime"
done

touch "$runtime/.ready"
