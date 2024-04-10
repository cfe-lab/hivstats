#! /bin/sh

set -xe

make reanalyze
test -z "$(git status --porcelain)"
