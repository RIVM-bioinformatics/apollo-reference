#!/bin/bash
# convert the bam of a Whole Genome Alignment of two assemblies into vcf format
#
# requires samtools (installed in the environment)
# requires k8 (and downloads it when not there)
# requires minimap2/misc/paftools.js (and downloads it when not there)

help_text="$(basename $0) <BAM> <FASTA>; writes uncompressed VCF to stdout"
if [ $(echo "$@" | grep -wcP "(\-h|\-\-help)") -eq 1 ]; then
  echo "$help_text"
  exit 0
elif [ $# -eq 2 ] && [ -f $1 ] && [ -f $2 ]; then
  bam=$(readlink -f $1)
  fasta=$(readlink -f $2)
  paf=/tmp/$(basename $bam).sorted.paf
else
  echo "$help_text"
  exit 0
fi

cd $(dirname $(readlink -f $0))
mkdir -p utils

if command -v paftools.js; then
  exe_paftools=paftools.js
elif [ ! -f utils/paftools.js ]; then
  # requires paftools.js
  wget https://raw.githubusercontent.com/lh3/minimap2/v2.30/misc/paftools.js -P utils
else
  exe_paftools=utils/paftools.js
fi


if command -v k8; then
  exe_k8=k8
elif [ ! -f utils/k8 ]; then
  # requires k8 (assuming linux)
  cd utils
  version=1.2
  wget -q https://github.com/attractivechaos/k8/releases/download/v$version/k8-$version.tar.bz2
  tar -jxvf k8-$version.tar.bz2
  cp k8-$version/k8-x86_64-Linux k8
  chmod a+x k8
  cd ../
  exe_k8=utils/k8
else
  exe_k8=utils/k8
fi

# BAM -> PAF -> VCF
samtools view -h --remove-flags 0x900 $bam | $exe_k8 $exe_paftools sam2paf - | sort -k6,6 -k8,8n > $paf
$exe_k8 $exe_paftools call -f $fasta $paf
exit 0
