#!/usr/bin/env bash
# Publish the built dataset to the main branch as a new, parentless commit.
#
#   ./publish.sh v1.0.0
#
# This repo keeps two unrelated histories in one place:
#
#   scripts  the generators (this branch). Normal linear history.
#   main     the data, and nothing else. One parentless commit per release.
#
# Each release replaces main rather than adding to it, so versions never
# compound: old ones survive only through their tags and can be pruned by
# deleting those tags, with no rewriting. See README.md.
#
# Nothing is pushed. The commands to do that are printed at the end.
set -euo pipefail
cd "$(dirname "$0")"

TAG=${1:-}
WT=.publish

die() { echo "error: $*" >&2; exit 1; }

[[ -n "$TAG" ]] || die "usage: ./publish.sh <tag>   e.g. ./publish.sh v1.0.0"
[[ "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || die "tag must look like v1.2.3, got '$TAG'"

[[ "$(git rev-parse --abbrev-ref HEAD)" == "scripts" ]] \
  || die "publish from the scripts branch (git switch scripts)"

# The generator commit is recorded in MANIFEST.md, so it has to mean something
[[ -z "$(git status --porcelain)" ]] \
  || die "commit or stash your changes first -- MANIFEST.md records the generator commit"

git rev-parse -q --verify "refs/tags/$TAG" >/dev/null \
  && die "tag $TAG already exists; bump the version"

[[ -d data ]] || die "no data/ -- run ./build.sh first"

# Samplesheets have the tag baked into every URL, so the two must agree or the
# published sheets will point at a version that does not exist.
CONFIG_TAG=$(python3 -c 'import sys; sys.path.insert(0, "scripts"); import config; print(config.TAG)')
REPO=$(python3 -c 'import sys; sys.path.insert(0, "scripts"); import config; print(config.REPO)')
[[ "$CONFIG_TAG" == "$TAG" ]] \
  || die "scripts/config.py has TAG='$CONFIG_TAG' but you asked to publish '$TAG'.
       Set TAG in scripts/config.py, run ./build.sh sheets manifest, commit, retry."

echo "==> validating"
./build.sh validate >/dev/null || die "validation failed; not publishing"
./build.sh manifest >/dev/null

echo "==> staging $TAG into a fresh root commit on main"
git worktree remove --force "$WT" 2>/dev/null || true
rm -rf "$WT"
git worktree prune
git worktree add --quiet --detach "$WT"
git -C "$WT" checkout -q --orphan _release
git -C "$WT" rm -rq --cached . >/dev/null 2>&1 || true
find "$WT" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +

rsync -a --exclude='.*' data "$WT"/
cp MANIFEST.md "$WT"/MANIFEST.md
python3 - "$WT/README.md" "$REPO" "$TAG" <<'PY'
import sys
from pathlib import Path
out, repo, tag = sys.argv[1], sys.argv[2], sys.argv[3]
text = Path("assets/README.data.md").read_text()
Path(out).write_text(text.replace("OWNER/REPO", repo).replace("/TAG", f"/{tag}"))
PY

git -C "$WT" add -A
git -C "$WT" commit -qm "Dataset $TAG

Built from generator commit $(git rev-parse HEAD) on branch scripts."
SHA=$(git -C "$WT" rev-parse HEAD)

git branch -f main "$SHA"
git tag -a "$TAG" -m "Dataset $TAG" "$SHA"

git -C "$WT" checkout -q --detach "$SHA"
git branch -D _release >/dev/null
git worktree remove --force "$WT"
git worktree prune

FILES=$(git ls-tree -r --name-only main | wc -l | tr -d ' ')
SIZE=$(git cat-file --batch-all-objects --batch-check='%(objectsize:disk)' 2>/dev/null \
       | awk '{t+=$1} END {printf "%.0f", t/1048576}')

cat <<EOF

published locally:
  main   -> $(git rev-parse --short main)  ($FILES files, parentless)
  tag    -> $TAG
  repo   -> $REPO
  .git   -> ${SIZE} MB

push both:
  git push --force origin main
  git push origin $TAG

then verify a file is actually served:
  curl -sI https://raw.githubusercontent.com/$REPO/$TAG/data/references/ecoli/genome.fasta | head -1
EOF
