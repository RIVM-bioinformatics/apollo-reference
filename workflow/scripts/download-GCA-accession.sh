#!/bin/bash
# download a NCBI WGS accession and the metadata we need

. $(dirname $(readlink -f $0))/build-helpers.sh
set +eu
accession=$1 || true
outdir=$2 || true # literal output directory (a subdirectory in  refdir ...)
set -eu

if [ $(echo $accession | grep -cPo "^GC[AF]_\d{9}\.\d+$") -eq 0 ]; then
  echo "exit: not a GC[AF]_\d{9}\.\d+ accession"
  exit 1;
elif [ -z "$outdir" ] || [ ! -d $outdir ]; then
  echo "exit: required existing \$outdir as 2th argument"
  exit 1;
fi
validate_refdir $(dirname $outdir)

# not used (yet) force re-downloading existing data
forced=0

## Requires (online) eutils suite;
## Eutils can be downloaded as stand-alone tools to, but for our use-case this suits the needs
#eutils=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
#efetch=$eutils/efetch.fcgi
#esearch=$eutils/esearch.fcgi
#datasets=$(dirname $(readlink -f $0))/../../utils/datasets
#samtools=samtools

# !important! (incrementally) create files in output directory
cd $outdir

if [ ! -f $accession.md5sum.txt ] || [ $forced -eq 1 ]; then
  # download
  echo "# start downloading: $accession ..."
  $datasets download genome accession $accession --filename $accession.zip
  # extract and copy the required files
  unzip -o $accession.zip
  cp ncbi_dataset/data/assembly_data_report.jsonl $accession.assembly_data_report.jsonl
  cp md5sum.txt $accession.md5sum.txt
  cp $(find ncbi_dataset/data/$accession -name "$accession*.fna") .
  # cleanup
  rm -f README.md
  rm -f md5sum.txt
  rm -rf ncbi_dataset
  rm -f $accession.zip
  # now obtain taxId and download its record
  # view human-readable
  # jq --color-output . $accession.assembly_data_report.jsonl
  json=$accession.assembly_data_report.jsonl
  taxId=$(cat $json | python3 -c "import sys, json; print(json.load(sys.stdin)['organism']['taxId'])")
  curl -s $efetch?db=taxonomy\&id=$taxId\&retmode=xml > $accession.$taxId.xml
fi
#ls -al $accession.md5sum.jsonl
ls -al $accession*

echo "# EOF=1 [$(basename $0)]"
exit 0








