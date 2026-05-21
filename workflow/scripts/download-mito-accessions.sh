#!/bin/bash
# download a/all NCBI nuccore corresponding mitochondrial accessions
# - includes prioritization to "*the* mitochondrion in case 1:many relationship
# - reports selected.mitochondria.tsv stating this prioritization
# - reports included.mitochondria.tsv stating co-deposited mitochondria in the genomic WGS entries

. $(dirname $(readlink -f $0))/build-helpers.sh
set +eu
refdir=$1 || true
validate_refdir $refdir
set -eu

# not used (yet) force re-downloading / recreation of existing data
recreate=0
cd $refdir

function report_tsv() {
  # DRY helper function to rewrite systematically formatted line to selected.mitochondria.tsv
  # Mind this same function is (mis)used for provided.mitochondria.tsv too;
  # the generated file is mv'ed lateren
  sed 's|NUCCORE/||g' \
  | sed 's/.fai:/\t/' | cut -f 1-4 | sed 's/^\(.\+\)\(.mito..\+\)$/\1\t\0/' | cut -f 2- \
  | tee /dev/stderr >> selected.mitochondria.tsv
}

function select_the_mitochondrion() {
  accession=$1
  any_available=$(ls NUCCORE/$accession.mito.*.fa 2>/dev/null | wc -l)
  if [ $any_available -ge 3 ]; then
    # prioritize to *THE* mitochondrion
    ls -altr NUCCORE/$accession.mito.*.fa | cat -n | tail -n 5
    # delete entries with >1 accession
    ls NUCCORE/$accession.mito.*.fa.fai \
      | xargs -i wc -l {} | awk '{ if ($1==1) { print $2 } }' \
      | xargs -i grep -HP "\d" {} \
      | awk '{a[NR]=$0; v[NR]=$2} END {m=int((NR+1)/2); print NR"\t"a[m]}' \
      | report_tsv
  elif [ $any_available -eq 2 ]; then
    # pick the largest one
    grep -HP "\d" NUCCORE/$accession.mito.*.fa.fai  | sort -gr -k 2 | head -n 1 | awk '{ print 2"\t"$0 }' \
      | report_tsv
  elif [ $any_available -eq 1 ]; then
    # there's just one
    grep -HP "\d" NUCCORE/$accession.mito.*.fa.fai | awk '{ print 1"\t"$0 }' \
      | report_tsv
  fi
  }

# re-generated the output files:
# - selected.mitochondria.tsv
# - included.mitochondria.tsv
# - provided.mitochondria.tsv     <-- TODO!
rm -f included.mitochondria.tsv
rm -f selected.mitochondria.tsv
rm -f provided.mitochondria.tsv
touch selected.mitochondria.tsv
touch provided.mitochondria.tsv

# first, download provided MT accessions
providedtxt=provided-MT-accessions.txt
if [ -f $providedtxt ]; then
  sed 1d $providedtxt | while read pair; do
    IFS=, read gca mt <<< "$pair"
    echo $gca $mt
    fa=NUCCORE/$gca.mito.$mt.fa
    if [ ! -f $fa ]; then
      curl -s $efetch?db=nuccore\&id=${mt}\&rettype=fasta > $fa
      samtools faidx $fa
      ls -al $fa
      head -n 2 $fa
      sleep 2
    fi
    grep -HP "\d" $fa.fai | awk '{ print 1"\t"$0 }' \
    | report_tsv
  done
  # !important! `report_tsv` writes to selected.mitochondria.tsv,
  # so we need to mv this file (and re-touch if for later appending)
  mv selected.mitochondria.tsv provided.mitochondria.tsv
  rm -f $providedtxt
fi

# loop over all the accessions and (try to) obtain corresponding mitochondria
find WGS -name "GC*.*.xml" | while read -r fname; do
  taxId=$(echo $fname | grep -Po "\.\d+\.xml$" | cut -f 2 -d'.')
  accession=$(echo $fname | awk -F'/' '{ print $NF }' | sed 's/\.[0-9]\+\.xml$//')
  echo "#" $accession $taxId
  min_nt_size=15000
  max_nt_size=100000
  # !important! realize capped at RetMax=20
  curl -s $esearch?db=nuccore\&term=txid"$taxId"+AND+mitochondrion%5Bfilter%5D+AND+"$min_nt_size"%3A"$max_nt_size"%5BSLEN%5D \
    | grep -Po "<Id>\d+</Id>" | grep -Po "\d+" | while read -r recId;
  do
    fa=NUCCORE/$accession.mito.$recId.fa
    if [ ! -f $fa ]; then
      curl -s $efetch?db=nuccore\&id=$recId\&rettype=fasta > $fa
      samtools faidx $fa
      ls -al $fa
      head -n 2 $fa
      sleep 2
    fi
  done
  # delete entries with >1 accession
  ls NUCCORE/$accession.mito.*.fa.fai 2>/dev/null \
    | xargs -i wc -l {} | awk '{ if ($1>=2) { print $0 } }' \
    | awk '{ print $0" deleting since >1 accession" }' | tee /dev/stderr \
    | awk '{ print $2 }' | xargs -i rm {}
  # select *THE* mitochondrion
  select_the_mitochondrion $accession
done

echo "# get information on mitochondria co-deposited in the assemblies"
# !important! genomic *.fai files needed
ls WGS/GC*_*.fna | xargs -i samtools faidx {}

find WGS -name "GC*.*.xml" | while read -r fname; do
  accession=$(echo $fname | awk -F'/' '{ print $NF }' | sed 's/\.[0-9]\+\.xml$//')
  assembly=$(ls WGS/$accession*.fna)
  num=$(grep "^>" $assembly | grep -cwiP "(mitochondrion|mitochondrial)" || true)
  grep "^>" $assembly | grep -wiP "(mitochondrion|mitochondrial)" \
    | tee /dev/stderr \
    | cut -f 2 -d'>' | cut -f 1 -d' ' | xargs -i grep -wH "^"{} $assembly.fai \
    | sed 's/.fai:/\t/' | cut -f 1-3 | awk '{ print "'$accession'\t'$num'\t"$0 }' \
    | sed 's|\tWGS/|\t|'
done > included.mitochondria.tsv

echo "# summary of all downloaded mitochondrial files:"
ls -altr NUCCORE/*.mito.*.fa 2>/dev/null | cat -n | tail -n 5
echo "# selected mitochondria files:"
cat selected.mitochondria.tsv
echo "# provided mitochondria (files):"
cat provided.mitochondria.tsv
echo "# included mitochondria (files):"
cat included.mitochondria.tsv
echo "# EOF=1 [$(basename $0)]"
exit 0

