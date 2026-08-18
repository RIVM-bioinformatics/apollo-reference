#!/bin/bash
# update (by downloading all reference assemblies listed in the provided TSV to the output directory
#
# This pipeline handles all what's shown in files/apollo-reference.drawio(.png)
#

help_text="$(basename $0) --out [/PATH/TO/OUT/DIR] [ --download-fastq ]"

# TODO: once hard-coding to Apollo_clinical_Candida_species.xlsx is released, update help_text and parse_args() 
#help_text="$(basename $0) --tsv [/PATH/TO/XLSX2CSV/APPOINTED/references.tsv] --out [/PATH/TO/OUT/DIR]"
#tsv=""

# arguments obtained from parse_args()
outdir=""
verbose=false
download_fastq=false
download_fastq_flag=""

# variables stating (sub)folders in this repo 
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
      --download-fastq)
        download_fastq=true
        download_fastq_flag="--download --recreate"
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


function logf() {
  # based on verbose, silence the output of some of the scripts
  if [ $verbose == true ]; then
    cat
  else
    tail -n 1
  fi
}

function summarize_downloaded_files() {
  # summarize downloaded files;
  # use quoted argument as input e.g. "$outdir/WGS/*.fna"
  filepattern=$1
  ls -tr $filepattern | xargs -i ls -al --time-style=+%Y%m%d  {} \
    | cut -f 3- -d' ' | cat -n | tail -n 5
  # log empty files / thus erroneous
  find $filepattern -type f -size 0 \
    | awk '{ print "# ERROR: zero-bytes file "$1 }' 1>&2
  }


function convert_xlsx_to_tsv() {
  # convert the appointed xlsx into appointed tsv file

  function get_accession_md5sum() {
    # helper function to termine unique accession fingerprint of tsv file (from xlsx)
    local tsvfile=$1
    csvcut -t -c "$colname_ref_accession","$colname_MT_accession" $tsvfile | grep GCA | md5sum | awk '{ print $1 }'
    }

  function regenerate_tsv_from_xlsx() {
    xlsx2csv $xlsx -d '\t' > $tsv
    # fix empty rows/cols: typical by-product of xlsx2csv i.c.w. Excel/LibreOffice combination
    # remove 100% empty columns
    while [ $(awk -F'\t' '{ print $NF }' $tsv  | sed '/^$/d' | wc -l) -eq 0 ]; do
      sed 's/\t$//' $tsv > $tsv.tmp
       mv $tsv.tmp $tsv
    done
    wc -l $tsv
    # remove 100% empty rows
    cat $tsv | sed '/^\s*$/d' > $tsv.tmp
    mv $tsv.tmp $tsv
    wc -l $tsv
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

# mind duplicated variable name in build-helpers.sh
refdata_tsv=$outdir/reference_assembly_data.tsv

echo "# xlsx          : $xlsx"
echo "# (input_)tsv   : $tsv"
echo "# refdata_tsv   : $refdata_tsv"
echo "# download_fastq: $download_fastq"

# create output subdirectories
mkdir -p $outdir/WGS
mkdir -p $outdir/SRR
mkdir -p $outdir/NUCCORE
mkdir -p $outdir/mmidx
mkdir -p $outdir/refs
mkdir -p $outdir/sam

# 1. convert XSLX to TSV (and pretty-show it)
convert_xlsx_to_tsv;
head $tsv | csvlook -t -I --snifflimit 0

# !important! need to keep the tsv in the outdir too (for build_helpers.sh input_tsv)
cp $tsv $outdir/$(basename $tsv)

# TODO: add validation of the xlsx with the available python code stated in this repo ....

# 2a. convert into plain list of accessions(.txt). This is the "queue" file for this workflow
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
  $scriptsdir/download-GCA-accession.sh $accession $outdir/WGS
done | logf

# 3b. list downloaded assembly files
summarize_downloaded_files "$outdir/WGS/*.fna"
echo "# EOF=1 [download-GCA-accession.sh]"

# 4. download mitochondria (given the TaxID of the assemblies
$scriptsdir/download-mito-accessions.sh $outdir | grep "^#"

# 4b. list downloaded mitochondrion files
summarize_downloaded_files "$outdir/NUCCORE/*.fa"
echo "# EOF=1 [download-mito-accessions.sh]"

### 5a. link accessions to SRR accessions (to generate self-test-data)
##$scriptsdir/link-accession-to-SRR.sh $outdir

# 5b. and subsequently download these SRR PE fastq datasets
$scriptsdir/link-accession-to-SRR.sh $outdir $download_fastq_flag | logf

# 6. generate "final" reference assembly data sheet
python3 $scriptsdir/generate_reference_assembly_dataframe.py $outdir | tee /dev/stderr > $refdata_tsv
# make sure the refdata_tsv - corresponding to the xlsx - is copied into the repo itself too!
cp $refdata_tsv $datadir/$(basename $refdata_tsv)
ls -al $refdata_tsv $datadir/$(basename $refdata_tsv)

# 7. generate per-species reference (including mitochondrion, if applicable)
$scriptsdir/build-species-references.sh $outdir

summarize_downloaded_files "$outdir/refs/*.fa"
echo "# EOF=1 [build-species-references.sh]"

# 8. generate (in case any update) new index for species identification
echo "# Realize there can/will be warnings like these:"
echo "# warning: MT accession xxxxxxx.y already included in other clade/WGS"
echo "# This means an MT accession is used twice for (multi-clade) references."
echo "# This needs to get corrected for in the minimap index"
$scriptsdir/build-identify-species-index.sh $outdir

echo "# EOF=1 [$(basename $0)]"
exit 0


}