#! /bin/sh

export DEBIAN_FRONTEND=noniteractive

set -xe

apt-get update
apt-get install -y mafft ncbi-blast+
apt-get install -y python3 python3-pip
apt-get install -y git make curl tar

curl -LsSf https://astral.sh/uv/install.sh > "/tmp/uv-install.sh"
sh "/tmp/uv-install.sh"
cp -v -f -- ~/.local/bin/uv ~/.local/bin/uvx /bin

uv sync
