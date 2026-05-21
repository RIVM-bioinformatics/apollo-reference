#!/bin/bash
# build genomic+MT DNA references for all species/accessions registered in the provided .tsv file

. $(dirname $(readlink -f $0))/build-helpers.sh
# target directory definition & validation
set +eu
refdir=$1 || true
validate_refdir $refdir
validate_input_tsv $input_tsv;
outdir=$refdir/refs
set -eu

# ------------------------------------------------------------------------------------------------- #
# TODO: implement use-case >1 MT_accession in WGS file. Not implemented yet (neither observed yet)
#       proposed solution: write a $outfa.yml file with this format:
#       ---
#       mitochondrial_accession:
#         - accession1
#         - accession2
#         - ...
# ------------------------------------------------------------------------------------------------- #
# TODO: implement check if there are non-annotated (as such) MT accessions in a genome
#       proposed solution: (blastn) provided/external MT to WGS
#                          in case hits found: .... nasty ....
#
#       Here, we could decide that the choice for reference genome +/- mitochondria is
#       the responsability of the biological researcher, not of the pipeline maintainer
#
# ------------------------------------------------------------------------------------------------- #

# not implemented as argument: toggle to 'true' to force re-creation of existing data
# !important! Realize end of this scripts writes-protects all created references
#             So, recreate=true will fail yelling:
#             cp: cannot create regular file '<fasta>': Permission denied
#             For the moment, I consider this acceptable behaviour, since overwriting
#             existing files is not an existing use-case.
#             If ever fasta editing (e.g. changing case, manipulating fasta headers, ...)
#             will be required for downstream tools acting on (indices of) these references,
#             simply toggle recreate=true and manually (once) set chmod a+w /path/to/.../refs/*.fa*
recreate=false

function generate_fasta_indices() {
  # generate .fai and .fal (for bedtools -g) indices for generated fasta
  fasta=$1
  samtools faidx $fasta
  cut -f 1-2 $fasta.fai > $fasta.fal
  ## bed-index of fasta isn't needed (yet?)
  #cut -f 1-2 $fasta.fai | awk '{ print $1"\t0\t"$2 }' > $fasta.bed
}


declare -a reference_array=()
convert_tsv_to_nested_reference_array $input_tsv reference_array

for ((i=0; i<${#reference_array[@]}; i+=5)); do
  read -r WGS_accession MT_accession MT_assembly MT_source ref_fasta <<< "${reference_array[@]:i:5}"
  wgs=$(find $refdir/WGS -name "$WGS_accession*_genomic.fna")

  if [ ! -f $ref_fasta ] || [ $recreate == true ]; then
    if [ "$MT_source" == "provided" ]; then
      # manually provided mitochondrion accession: concatenate WGS + mitochondrial fasta
      echo "generating [$MT_source] --> $ref_fasta"
      cat $wgs $refdir/NUCCORE/$MT_assembly > $ref_fasta
      generate_fasta_indices $ref_fasta
      # generate associated *.mitochondrion.fal index
      cut -f 1-2 $refdir/NUCCORE/$MT_assembly.fai > $ref_fasta.mitochondrion.fal
    elif [ "$MT_source" == "external" ]; then
      # external mitochondrion: concatenate WGS + mitochondrial fasta
      echo "generating [$MT_source] --> $ref_fasta"
      cat $wgs $refdir/NUCCORE/$MT_assembly > $ref_fasta
      generate_fasta_indices $ref_fasta
      # generate associated *.mitochondrion.fal index
      cut -f 1-2 $refdir/NUCCORE/$MT_assembly.fai > $ref_fasta.mitochondrion.fal
    elif [ "$MT_source" == "included" ]; then
      # internally defined mitochondrion: copy WGS fasta
      echo "generating [$MT_source] --> $ref_fasta"
      cp $wgs $ref_fasta
      generate_fasta_indices $ref_fasta
      # generate associated *.mitochondrion.fal index
      grep "^>" $ref_fasta | grep -wiP "(mitochondrion|mitochondrial)" \
        | awk '{ print $1 }' | tr -d ">" \
        | xargs -i grep -m 1 -w "^"{} $ref_fasta.fal \
        > $ref_fasta.mitochondrion.fal
    elif [ "$MT_source" == "NA" ] || [ "$MT_source" == "ND" ]; then
      # be forward-compatible in case NA/ND difference for defined MT accession is introduced
      # no mitochondrion defined neither available in NUCCORE: copy WGS fasta
      echo "generating [ND] --> $ref_fasta"
      cp $wgs $ref_fasta
      generate_fasta_indices $ref_fasta
    else
      echo "error: NotImplementedError(MT_source=$MT_source)"
      exit 1
    fi
  fi
done

# make all files +r / -w for everybody
chmod a+r $outdir/*
chmod a-w $outdir/*
echo "# summary of all generated reference files:"
ls -altr $outdir/*.fa 2>/dev/null | cat -n | tail -n 5
echo "# summary of origin of mitochondrial data:"
ls -altr $outdir/*.fa 2>/dev/null | awk -F'__' '{ print $2 }' | sort | uniq -c
echo "# EOF=1 [$(basename $0)]"
exit 0

