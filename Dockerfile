
FROM ubuntu:22.04

WORKDIR /workspaces/hivstats
COPY . .

RUN sh src/install-dependencies.sh

ENTRYPOINT make reanalyze
