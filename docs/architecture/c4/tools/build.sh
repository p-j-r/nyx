#!/usr/bin/env bash
# Rebuilds the whole C4 set: reads the Rust source and spec/components.spec, writes the Structurizr fragments and the code pages,
# validates the workspace with the Structurizr CLI, exports Mermaid and PlantUML, and writes the pages and the viewer.
#
#   bash tools/build.sh                       # in the repository: docs/architecture/c4 next to docs/architecture/workspace.dsl
#   bash tools/build.sh --standalone          # in the standalone folder: structurizr/workspace.dsl, structurizr/adr
#
# Needs: perl, Java, and the Structurizr CLI release (https://github.com/structurizr/cli/releases, unzipped).
# Point STRUCTURIZR_CLI at its folder (the one with lib/); default ../../../structurizr-cli/unpacked relative to the repository.
# CRATES points at the repository's crates/ folder (default: three folders up from here for the repository layout).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mode=repo; [ "${1:-}" = "--standalone" ] && mode=standalone
if [ "$mode" = repo ]; then
  dsl="$here/../workspace.dsl"; model="model"; adr="../adr"; crates="${CRATES:-$here/../../../crates}"
else
  dsl="$here/structurizr/workspace.dsl"; model="structurizr/model"; adr="structurizr/adr"; crates="${CRATES:?set CRATES to the repository crates folder}"
fi
cli="${STRUCTURIZR_CLI:-}"
if [ -z "$cli" ]; then for c in "$here/../../../../structurizr-cli/unpacked" "$here/../structurizr-cli/unpacked"; do [ -d "$c/lib" ] && cli="$c" && break; done; fi
[ -n "$cli" ] && cli="$(cd "$cli" && (pwd -W 2>/dev/null || pwd))"
[ -d "$cli/lib" ] || { echo "Structurizr CLI not found at $cli (set STRUCTURIZR_CLI)"; exit 2; }
sz() { java -cp "$cli/lib/*" com.structurizr.cli.StructurizrCliApplication "$@"; }

cd "$here"
perl tools/gen.pl --root . --crates "$crates" --model "$model" --adr "$adr" --dsl "$dsl"
sz validate -workspace "$dsl"
rm -rf export; mkdir -p export
sz export -workspace "$dsl" -format mermaid  -output export/mermaid  > /dev/null
sz export -workspace "$dsl" -format plantuml -output export/plantuml > /dev/null
sz export -workspace "$dsl" -format json     -output export/json     > /dev/null
perl tools/gen.pl --pages --root . --crates "$crates" --model "$model" --adr "$adr" --dsl "$dsl"
rm -f reports/.codediagrams.tmp
echo "built: $(ls export/mermaid | wc -l) views exported"
