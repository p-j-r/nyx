#!/usr/bin/env bash
# Fails when a tracked file matches no row of docs/architecture/CODE_TO_ADR.md, or when a row names an ADR that does not exist.
# Run from anywhere:  bash docs/architecture/check_code_to_adr.sh
set -u
cd "$(dirname "$0")/../.." || exit 2
map=docs/architecture/CODE_TO_ADR.md
fail=0
set -f   # brace patterns are expanded by eval below, but never against the file system

rows=$(sed -n '/<!-- map:begin -->/,/<!-- map:end -->/p' "$map" | grep '^| `' | sed -E 's/^\| `([^`]*)` \| ([^|]*) \|.*/\1\t\2/')

# patterns (brace-expanded once) and every ADR number named in the table
patterns=()
adrs=""
while IFS=$'\t' read -r glob nums; do
  [ -z "$glob" ] && continue
  for g in $(eval "echo $glob"); do patterns+=("$g"); done
  adrs="$adrs,$nums"
done <<< "$rows"

for n in $(printf '%s\n' "$adrs" | tr ',' '\n' | tr -d ' ' | grep -E '^[0-9]+$' | sort -un); do
  f=$(printf 'docs/architecture/adr/%04d-*.md' "$n")
  set +f; ls $f >/dev/null 2>&1 || { echo "CODE_TO_ADR names ADR $n but $f does not exist"; fail=1; }; set -f
done

if command -v jj >/dev/null 2>&1 && jj file list >/dev/null 2>&1; then
  files=$(jj file list | tr '\\' '/')
elif git ls-files >/dev/null 2>&1 && [ -n "$(git ls-files | head -1)" ]; then
  files=$(git ls-files)
else
  files=$(find . -type f -not -name '.git' -not -path './.git/*' -not -path './.jj/*' -not -path '*/target/*' -not -path './dist/*' -not -path './runs/*' | sed 's|^\./||')
fi

missing=0
total=0
while IFS= read -r f; do
  [ -z "$f" ] && continue
  total=$((total+1))
  hit=0
  for g in "${patterns[@]}"; do
    # shellcheck disable=SC2053
    if [[ $f == $g ]]; then hit=1; break; fi
  done
  if [ $hit = 0 ]; then echo "no row for: $f"; missing=$((missing+1)); fi
done <<< "$files"

[ $missing -gt 0 ] && fail=1
echo "$total files checked, $missing without a row, ${#patterns[@]} patterns"
exit $fail
