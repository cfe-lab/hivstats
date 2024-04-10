#! /bin/sh

set -xe

docker build -t hivstats .
rm -rf outputs
docker run -v "$PWD":"/workspaces/hivstats" hivstats
test -z "$(git status --porcelain)"
