#!/bin/bash
# build minimap2 index corresponding to all references registered in the provided ***references.tsv file
#
# !important! if adjusting, realize exterior apollo_reference expects this script to report
#             as FINAL unquoted line the full path to the correspoding minimap2 assembly
#             By calling this script, one can retrieve (optionally build first)
#             the corresponding species identification index given the $input_tsv
#
#
. $(dirname $(readlink -f $0))/build-helpers.sh
# target directory definition & validation
set +eu
refdir=$1 || true
validate_refdir $refdir
# overrule input_tsv to the provided refdir
input_tsv=$refdir/reference_assembly_data.tsv
validate_input_tsv $input_tsv;
set -eu

# generate the (unique) concatenated minimap index name based on accessions in $input_tsv
datasetUID=$(generate_reference_accession_dataset_UID $input_tsv)
mmidx=$refdir/mmidx/apollo-species-refs.$datasetUID.mmidx

if [ 1 -eq 0 ] && [ -f $mmidx ]; then
  echo "# expected minimap2 index [$datasetUID] already exists!"
else
  echo "# create minimap2 index from concatenated reference fasta files ..."
  declare -a reference_array=()
  convert_tsv_to_nested_reference_array $input_tsv reference_array
  mmidx_tmp_fasta=/tmp/concatenated.$datasetUID.fasta
  rm -f $mmidx_tmp_fasta
  touch $mmidx_tmp_fasta
  ls -al $mmidx_tmp_fasta
  # 1. Add all the genome sequences themselves
  #    Since all coming from NCBI/WGS, accession **should** be unique
  for ((i=0; i<${#reference_array[@]}; i+=5)); do
    read -r WGS_accession MT_accession MT_assembly MT_source ref_fasta <<< "${reference_array[@]:i:5}"
    wgs_fasta=$(find $refdir/WGS -name "$WGS_accession*_genomic.fna")
    cat $wgs_fasta >> $mmidx_tmp_fasta
  done
  ls -al $mmidx_tmp_fasta
  # 2. Add custom mitochondria
  for ((i=0; i<${#reference_array[@]}; i+=5)); do
    read -r WGS_accession MT_accession MT_assembly MT_source ref_fasta <<< "${reference_array[@]:i:5}"
    wgs_fasta=$(find $refdir/WGS -name "$WGS_accession*_genomic.fna")
    if [ "$MT_source" == "included" ]; then
      skip=1
    elif [ "$MT_assembly" != "NA" ]; then
      if [ $(grep -c "^>$MT_accession" $mmidx_tmp_fasta) -eq 0 ]; then
        mit_fasta=$refdir/NUCCORE/$MT_assembly
        cat $mit_fasta >> $mmidx_tmp_fasta
      else
        # case multi-clade reference all defaulting to the same mitochondrion
        echo "warning: MT accession $MT_accession already included in other clade/WGS" 1>&2
        continue
      fi
    fi
  done
  ls -al $mmidx_tmp_fasta
  # build index and remove fasta file (index alone will suffice)
  minimap2 -d $mmidx $mmidx_tmp_fasta
  # Build an accession,length TSV alike a *.fasta.fal (or the -g file used in bedtools);
  # realize (over)write-protection
  if [ ! -f $mmidx.fal ]; then
    minimap2 -a $mmidx 2>/dev/null | grep "^@SQ" | cut -f2,3 | sed 's/SN://g; s/LN://g' > $mmidx.fal
  fi
  # QC: all accessions in mmidx should be unique;
  # - prevents downstream htslib errors [W::sam_hdr_create], [E::sam_hrecs_update_hashes]
  # - prevents is_valid_headered_sam failing in apollo-match-reference
  if [ $(cut -f 1 $mmidx.fal | sort | uniq -c | awk '{ if ($1!=1) { print $0 } }' | tee /dev/stderr | wc -l) -ge 1 ]; then
    echo "error: duplicated (likely mitochondrion ...) accessions in mmidx [$mmidx]" 1>&2
    rm -f $mmidx_tmp_fasta 
    rm -f $mmidx
    rm -f $mmidx.fal
    exit 1
  fi
  rm -f $mmidx_tmp_fasta
fi

# settle modification rights;
# since apollo-reference will be future ISO-certified, existing (reference) files + mmidx can't be modified
find "$refdir" -type d -exec chmod a+rwx,a+t {} +
find "$refdir" -type f -exec chmod a+r-w {} +
find "$refdir" -maxdepth 1 -type f -exec chmod ug+w {} +

# this is the FINAL unquoted line reported by this script
ls $mmidx*
ls -al $mmidx* | awk '{ print "# "$0 }'
echo "# EOF=1 [$(basename $0)]"
exit 0


# -------------------------------------------------------------------------------------------
# blueprint of first POC and performance test
# -------------------------------------------------------------------------------------------
# refdir=/mnt/db/apollo/reference
# mmidx=$refdir/mmidx/apollo-species-refs.77d492a325d33dce8cb42e82a914ef1e.mmidx
# SRR_accession=SRR9007776
# fq1=$refdir/SRR/"$SRR_accession"_1.fastq.gz
# fq2=$refdir/SRR/"$SRR_accession"_2.fastq.gz
# time minimap2 -t 10 -ax sr $mmidx $fq1 $fq2 > alignment.sam
# roughy 15 minutes for 2.88M 151nt PE pairs, ~70x coverage (bit less, MT is much higher covered)
shared="samtools view  alignment.sam -F 2048 -F 256 -f 2 | cut -f 3 | sort | uniq -c | sort -g | tail -n 50"
paste <(eval $shared) \
      <(eval $shared | awk '{ print $2 }' | xargs -i grep -m 1 -w "^"{} data/*.fai | cut -f 1-2)
# based on this first test, seems to work pretty well!
# off-mappings occur at "hotspots":
# dirt-picked example being the highest observed wrong accession,
# being CM076439.1
# bedtools bamtobed -i alignment.sam | grep -w CM076439.1 | bedtools sort | bedtools merge -d -50 -c 3 -o count | awk '{ cov=$4/($3-$2); if (cov>0.15) { print $0"\t"($3-$2)"\t"cov}  }'
#
# Sortes-vs-unsorted SAM, in ViroConstrictor:
# #Florian:
# Hoeft inderdaad niet per-se. we werken eigenlijk exclusief met Bam files en een sorted bamfile is een stuk efficienter voor random-access met pysam.
# Theoretisch is het hier niet echt nodig puur op basis van wat we nu doen.
# Maar op deze manier doen we praktisch dezelfde operatie aan beide kanten wat makkelijker is om te onderhouden en
# bij de initiële implementatie wisten we nog niet zeker of we het alleen bij read-count per referentie lieten
# of dat we ook nog andere informatie uit de alignment files wouden halen waarvoor een sorted bam wel nuttig kon zijn

