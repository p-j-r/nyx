#!/usr/bin/env bash
# Builds the standalone copy of the C4 set in a folder of its own (default: ../../../../nyx-c4, a sibling of the repository folders),
# with its own Structurizr workspace, the ADRs it links to, the exports, the pages and the viewer. The folder does not need the
# repository except to regenerate (CRATES=<repository>/crates bash tools/build.sh --standalone).
#
#   bash tools/make_standalone.sh [destination]
#
# The destination is only ever added to or overwritten file by file; nothing in it is deleted (the export folder is rebuilt by
# build.sh, which removes its own export/ and nothing else). Remove an old copy yourself if you want a clean start.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
repo="$(cd "$here/../../.." && pwd)"
dest="${1:-$here/../../../../nyx-c4}"
mkdir -p "$dest"; dest="$(cd "$dest" && pwd)"
[ "$dest" = "$here" ] && { echo "destination is the source folder"; exit 2; }
mkdir -p "$dest/structurizr/model" "$dest/structurizr/adr" "$dest/spec" "$dest/tools"
cp "$here"/spec/components.spec "$dest/spec/"
cp "$here"/tools/gen.pl "$here"/tools/build.sh "$here"/tools/make_standalone.sh "$dest/tools/"
cp "$here"/MEMO.md "$dest/"
sed 's#(\.\./workspace\.dsl)#(structurizr/workspace.dsl)#; s#(\.\./adr/)#(structurizr/adr/)#' "$here/README.md" > "$dest/README.md"
cp "$repo"/docs/architecture/adr/*.md "$dest/structurizr/adr/"
sed 's#!include c4/model/#!include model/#' "$repo/docs/architecture/workspace.dsl" > "$dest/structurizr/workspace.dsl"
CRATES="$repo/crates" bash "$dest/tools/build.sh" --standalone
echo "standalone C4 set: $dest"
