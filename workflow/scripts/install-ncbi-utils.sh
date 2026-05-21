#!/bin/bash
# install NCBI tools needed in this package
# - datasets
# - fasterq-dump
# - eutils (esearch/efetch) suite is used via the web API.
cd $(dirname $(readlink -f $0))
mkdir -p utils

if [ -z "$(which utils/datasets)" ]; then
  curl -s https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets > utils/datasets
  chmod +x utils/datasets
fi
./utils/datasets --version

if [ -f ./utils/fasterq-dump ]; then
  version=3.3.0
  url=https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${version}/sratoolkit.${version}-ubuntu64.tar.gz
  progname=$(basename $url | sed 's/.tar.gz$//')
  wget -P utils $url
  tar -vxzf utils/$progname.tar.gz -C utils
  rm -f utils/$progname.tar.gz
  # temporarily solve with a symlink ...
  ln -s $(pwd)/utils/$progname/bin/fasterq-dump $(pwd)/utils/fasterq-dump
fi
./utils/fasterq-dump --version