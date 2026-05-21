#!/bin/bash
# update (by downloading all reference assemblies listed in the provided TSV to the output directory
#
# This pipeline handles all what's shown in files/apollo-reference.drawio(.png)
#

help_text="$(basename $0) --out [/PATH/TO/OUT/DIR]"

#help_text="$(basename $0) --tsv [/PATH/TO/XLSX2CSV/APPOINTED/references.tsv] --out [/PATH/TO/OUT/DIR]"
#tsv=""

outdir=""
verbose=false
thisdir=$(dirname $(readlink -f $0))
packagedir=$(dirname $thisdir)
datadir=$packagedir/data
scriptsdir=$thisdir/scripts

function parse_args() {
  # parse provided arguments to this script
  while [[ $# -gt 0 ]]; do
    case $1 in
      --verbose|-v)
        verbose=true
        shift
        ;;
      # TODO: needed once hardcoded xlsx is replaced by tsv input
      #--tsv)
      #  tsv="$2"
      #  shift 2
      #  ;;
      --out)
        outdir="$2"
        shift 2
        ;;
      --help|-h|*)
        echo "$help_text"
        exit 0
        ;;
    esac
  done

  # TODO: needed once hardcoded xlsx is replaced by tsv input
  #if [ -z "$tsv" ] || [ ! -f $tsv ]; then
  #  echo "error: --tsv /path/to/existing/file.tsv is required." >&2
  #  echo "$help_text"
  #  exit 1
  if [ -z "$outdir" ] || [ ! -d $outdir ]; then
    echo "error: --out /path/to/existing/outdir is required." >&2
    echo "$help_text"
    exit 1
  elif [ "$outdir" == "." ] || [ "$outdir" == "$thisdir" ]; then
    echo "error: cowardly refusing to let --out be this directory (don't pollute your repo!)" >&2
    echo "$help_text"
    exit 1
  elif [ "$outdir" == "$packagedir" ]; then
    echo "error: cowardly refusing to let --out be this repository (don't pollute your repo!)" >&2
    echo "$help_text"
    exit 1
  fi
}

# (naive) format definition in xlsx; get accessions from these columns
# In case non-compatible xlsx/sheet, workflow will crash at some point
# In apollo_reference/referencedata.py, python apollo_reference to validate the sheet/tsv is stated.
# So, once workflow is ported to run_pipeline.py, checks can be done pythonically
# For now, since (full) support is there only for 1 reference-dataset,
# files (and names) are in this repo and are hard-coded here.
xlsx=$datadir/Apollo_clinical_Candida_species.xlsx
tsv=$datadir/supported-reference-species.tsv
colname_ref_accession="Reference accession"
colname_MT_accession="Mitochondrion accession"

echo "# xlsx: $xlsx"
echo "# tsv : $tsv"

function convert_xlsx_to_tsv() {
  # convert the appointed xlsx into appointed tsv file

  function get_accession_md5sum() {
    # helper function to termine unique accession fingerprint of tsv file (from xlsx)
    local tsvfile=$1
    csvcut -t -c "$colname_ref_accession","$colname_MT_accession" $tsvfile | grep GCA | md5sum | awk '{ print $1 }'
    }

  function regenerate_tsv_from_xlsx() {
    xlsx2csv $xlsx -d '\t' > $tsv
    # remove 100% empty columns (by-product of xlsx2csv.py)
    while [ $(awk -F'\t' '{ print $NF }' $tsv  | sed '/^$/d' | wc -l) -eq 0 ]; do
      sed 's/\t$//' $tsv > $tsv.tmp
       mv $tsv.tmp $tsv
    done
    }

  if [ ! -f $xlsx ]; then
    echo "exit: expected $xlsx in this repo ..."
    exit 1
  elif [ -f $xlsx ] && [ -f $tsv ]; then
    xlsx2csv $xlsx -d '\t' > /tmp/$(basename $tsv).tmp
    md5sum_EXISTING=$(get_accession_md5sum /tmp/$(basename $tsv).tmp)
    md5sum_CURRENT=$(get_accession_md5sum $tsv)
    if [ "$md5sum_EXISTING" == "$md5sum_CURRENT" ]; then
      echo "# no updates in tsv [$tsv] versus xlsx file [$xlsx]"
      echo "# are you sure you want to continue?"
      echo "# sleep 10 ... [ terminate using <Ctrl+C> ]"
      sleep 10
    else
      regenerate_tsv_from_xlsx;
    fi
  elif [ -f $xlsx ] && [ ! -f $tsv ]; then
    regenerate_tsv_from_xlsx;
  fi
}

# parse arguments & start script
parse_args "$@"

# NOT NEEDED YET ...
#. $scriptsdir/build-helpers.sh

# create output subdirectories
mkdir -p $outdir/WGS
mkdir -p $outdir/SRR
mkdir -p $outdir/NUCCORE
mkdir -p $outdir/mmidx
mkdir -p $outdir/refs
mkdir -p $outdir/sam

# 1. convert XSLX to TSV (and pretty-show it)
convert_xlsx_to_tsv;
head $tsv | csvlook -t -I;
# !important! need to keep the tsv in the outdir too (for build_helpers.sh input_tsv)
cp $tsv $outdir/$(basename $tsv)

# 2a. convert into plain list of accessions(.txt)
awk -F'\t' '{ if ($3!=1 && $8!="RIVM") { print $0 } }' $tsv \
  | csvcut -t -c "$colname_ref_accession" | sed 1d \
  > $outdir/accessions.txt
wc -l $outdir/accessions.txt

# 2b. convert into lookup of provided-MT-accessions(.txt)
outtxt=$outdir/provided-MT-accessions.txt
awk -F'\t' '{ if ($3!=1 && $8!="RIVM") { print $0 } }' $tsv \
  | csvcut -t -c "$colname_ref_accession","$colname_MT_accession" \
  | sed '/,$/d' | sed 1d | awk 'BEGIN { print "WGS,MT" } { print $0 }' \
  > $outtxt
if [ $(wc -l < $outtxt) -eq 1 ]; then
  rm -f $outtxt
else
  wc -l $outtxt
fi

# input data ready now' move to output directory and start the download workflow
cd $outdir

# 3. download reference genome assemblies
cat accessions.txt | while read accession; do
  $scriptsdir/download-GCA-accession.sh $accession $outdir/WGS;
done

# 4. download mitochondria (given the TaxID of the assemblies
$scriptsdir/download-mito-accessions.sh $outdir

# 5a. link accessions to SRR accessions (to generate self-test-data)
$scriptsdir/link-accession-to-SRR.sh $outdir

# 5b. and subsequently download these SRR PE fastq datasets
$scriptsdir/link-accession-to-SRR.sh $outdir --download

if [ 1 -eq 0 ]; then
  # 6. generate "final" reference assembly data sheet
  python3 $scriptsdir/generate_reference_assembly_dataframe.py $outdir
fi

# 7. generate per-species reference (including mitochondrion, if applicable)
$scriptsdir/build-species-references.sh $outdir

# 8. generate (in case anu update) new identify species index
$scriptsdir/build-identify-species-index.sh $outdir

echo "# EOF=1 [$(basename $0)]"
exit 0

