#!/usr/bin/env bash
set -euo pipefail

[[ -f Korenev_ID.tex ]] || { echo "Korenev_ID.tex not found" >&2; exit 1; }

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

zip -j "$tmp_dir/Korenev_ID.zip" Korenev_ID.tex
mv -f "$tmp_dir/Korenev_ID.zip" Korenev_ID.zip