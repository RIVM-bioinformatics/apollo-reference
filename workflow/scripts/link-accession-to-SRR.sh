#!/bin/bash
# link species accessions to corresponding SRR PE Illumina fastq data
# - expects and only supports SRR PE Illumina fastq data

. $(dirname $(readlink -f $0))/build-helpers.sh
refdir=$1 || true
validate_refdir $refdir

# ---------------------------------------------------------------------------------- #
# TODO: map to SRR via SRS and download raw Illumina fastq (for benchmarking)
# ---------------------------------------------------------------------------------- #
# issue here breaking straightforward "just download the fastq"
# is that many of these projects are:
# - one2many: "umbrella" projects (multi-species WGS projects)
# - one2many: multiple SRR-runs for a single WGS-assembly
# - much data is PacBio, not simple PE Illumina
# Because of this, it is for now decided not to automate the download of this data any further
# than for the few for which the below snippet does work nicely

# generated output filename
accession2SRR=accession2SRR.tsv
recreate=$(echo "$@" | grep -c "\-recreate" || true)
download=$(echo "$@" | grep -c "\-download" || true)
cd $refdir

if [ $recreate -eq 1 ]; then
  wc -l $accession2SRR
  printf "accession\tSRR\n" > $accession2SRR
  cat accessions.txt | while read -r accession; do
    json=WGS/$accession.assembly_data_report.jsonl
    samnId=$(cat $json | python3 -c "import sys, json; print(json.load(sys.stdin)['assemblyInfo']['biosample']['accession'])" 2>/dev/null || true)
    samnId=${samnId:-NA}
    echo "#" $accession $samnId #1>&2
    if [ "$samnId" == "NA" ]; then continue; fi
    curl -s $esearch?db=sra\&term="$samnId"+AND+ILLUMINA%5BPLATFORM%5D+AND+PAIRED%5BLAYOUT%5D\&rettype=runinfo\&retmode=text \
      | grep -Po "<Id>\d+</Id>" | tee /dev/stderr | grep -Po "\d+" | head -n 1 \
      | xargs -i curl -s $efetch?db=sra\&id={}\&rettype=runinfo \
      | grep -Po "<Run>SRR\d+</Run>" | sed 's/<[^>]\+>//g' \
      | awk '{ print "'$accession'\t"$1 }' | tee /dev/stderr \
      >> $accession2SRR
    sleep 5
  done
fi
wc -l $accession2SRR

# download fastq from SRA
cat accessions.txt | while read -r accession; do
  srrId=$(grep -w $accession accession2SRR.tsv | cut -f 2)
  if [ -z "$srrId" ]; then
    continue
  elif [ -f SRR/$srrId"_1.fastq.gz" ]; then
    continue
  elif [ $download -eq 0 ]; then
    echo "skipping download of SRR/"$srrId"_1/2.fastq.gz; use $(basename $0) --download to enable"
    continue
  else
    echo "# start downloading $srrId"
    $fasterqdump $srrId --split-files --outdir SRR --temp /tmp
    echo "gzipping fastq ... patience please"
    #tar -czf SRR/$srrId"_1.fastq.gz" SRR/$srrId"_1.fastq" --remove-files &
    #tar -czf SRR/$srrId"_2.fastq.gz" SRR/$srrId"_2.fastq" --remove-files &
    gzip SRR/$srrId"_1.fastq" &
    gzip SRR/$srrId"_2.fastq" &
    wait
    ls -al SRR/$srrId*.fastq.gz
  fi
done

# log all files downloaded so far
ls -altr SRR/*.fastq.gz | cat -n | tail
echo "# EOF=1 [$(basename $0)]"
exit 0