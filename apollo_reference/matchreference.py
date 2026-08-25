#!/usr/bin/env python3

# Python Imports
import os
import re
import glob
import sys
import yaml
import shutil
from collections import Counter
from tabulate import tabulate
from pathlib import Path
from typing import List, Tuple, NamedTuple
import pandas as pd
import subprocess
import argparse

# Package Imports
try:
    from dataframes import read_reference_species_and_assembly_df
    from configuration import REFERENCEDATA_YAML as CONFIGFILE
except ModuleNotFoundError:
    from .dataframes import read_reference_species_and_assembly_df
    from .configuration import REFERENCEDATA_YAML as CONFIGFILE

# External package Imports (generic argparse helpers)
from rivm_ids_swc_argparseutils.actions import DynamicHelpTopicAction
from rivm_ids_swc_argparseutils.shortcuts import as_argparse_type
from rivm_ids_swc_argparseutils.formats.fasta import validate_fasta_file

# read from config where the downloaded sequence data + mmidx indices should reside
REPODIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DATA_DIR = yaml.safe_load(open(CONFIGFILE))['reference_data_dir']

# output file suffixes for the *references.yml YAML file
# !important!   apollo-mapping Snakefile implements a wildcard rule that
#               will break in a renaming to (eg) forced-references.yml"
#               Generally, any renaming here should be adjusted in apollo-mapping,
#               since this/these yamls are the target of a crucial checkpoint
SUFFIX_YAML_MATCHED = "-references.yml"
SUFFIX_YAML_FORCED = "-forced_references.yml"

# Match any/path/to/GCA_017309295.1.418086.xml
# This is the most general, generic and single file that contains
# the information (linked to a corresponding *.fai file with the accessions)
# to which reference and its corresponding taxonomy any contig belongs
REGEX_PATH_TAX_RECORD = re.compile(r"^(?P<pathprefix>.+/)?(?P<reference>[^/]+)\.(?P<taxid>\d+)\.xml$")

# ------------------------------------------------------------------------------------- #
# argument parser functions
# ------------------------------------------------------------------------------------- #

if not shutil.which("samtools"):
    msg = ("Error: samtools executable not found in PATH.\n"
        "       %s assumes this working-horse piece of software to be (always) installed.\n"
        "       Please install it via your package manager (e.g., 'conda install -c bioconda samtools' "
        "or 'apt install samtools') before running this application.") % os.path.basename(__file__)
    sys.exit(msg)

def naive_is_headered_sam_file(parser:argparse.ArgumentParser,arg:str) -> Path:
    """ Naively test if the provided argument is a headered SAM file """
    cmd = f"samtools view -H {arg}"
    if not os.path.isfile(arg):
        parser.error(f"not a file: {arg}")
    elif not os.path.splitext(arg)[-1] == ".sam":
        parser.error(f"not a SAM/BAM file: {arg}")
    elif subprocess.run(cmd, shell=True, capture_output=True, text=True).stderr:
        parser.error(f"not a (headered) SAM/BAM file: {arg}")
    elif not subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout:
        parser.error(f"not a headered SAM/BAM file: {arg}")
    return Path(arg)

def is_valid_writeable_location_prefix(parser:argparse.ArgumentParser,arg:str) -> str:
    """ """
    if os.path.isdir(arg):
        parser.error(f"prefix is a directory [{arg}]; expected /path/to/existing/dir/PREFIX")
    elif not os.path.isdir(os.path.dirname(arg)):
        parser.error(f"directory component of prefix [{arg}] is not an existing directory")
    # TODO: test if writeable
    return arg

def add_out_prefix_argument(parser:argparse.ArgumentParser):
    """ adds the required out_prefix argument to the parser """
    parser.add_argument(
        "out_prefix",
        type=lambda x: is_valid_writeable_location_prefix(parser, x),
        metavar="PREFIX",
        help="file path prefix for yaml and the various dataframes with results",
    )

def get_parser() -> argparse.ArgumentParser:
    """ Construct the ArgumentParser for the/this match-ref.py script """
    descr = "process a SAM file (minimapped PE fastq) in search for the best corresponding species (taxID) and assembly"
    parser = argparse.ArgumentParser(
        description=descr,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "sam_file",
        type=lambda x: naive_is_headered_sam_file(parser, x),
        metavar="SAM",
        help="(unsorted) headered SAM file",
    )
    add_out_prefix_argument(parser)

    _group = parser.add_argument_group('Arguments to exclude (mitochondrial) accessions that occur in the index')
    group = _group.add_mutually_exclusive_group()

    def split_comma_separated_list(value:str) -> List[str]:
        return value.split(',')

    group.add_argument(
        "--excl-mito-accessions","--excl-accessions",
        type=split_comma_separated_list,
        metavar="A,B[,..]",
        help="exclude comma-separated list of (mitochondrial) accessions",
    )

    group.add_argument(
        "--excl-mito-tsv",
        #type=lambda x: naive_is_headered_sam_file(parser, x),
        metavar="TSV",
        help="exclude all mitochondrial accessions listed in provided reference_assembly_data.tsv",
    )

    class DynamicHelpUsageAction(DynamicHelpTopicAction):
        def generate_content(self) -> str:
            help_text = open(os.path.join(REPODIR, "docs", "README-match-ref.txt")).read()
            return "\n%s\n" % (help_text % {'prog':parser.prog})

    parser.add_argument("--help-usage",
        action=DynamicHelpUsageAction,
        help=f"show detailed explanation of how to use (or bypass!) %(prog)s"
    )
    return parser


def get_forced_parser() -> argparse.ArgumentParser:
    """ Construct the ArgumentParser for the "hook" of bypassing match-ref.py (which generates a "forced" yaml) """
    descr = "Trigger the generation of <OUTPREFIX>-forced.yml file (which bypasses %(prog)s in a pipeline)"
    parser = argparse.ArgumentParser(description=descr)

    add_out_prefix_argument(parser)

    _group = parser.add_argument_group('provide reference assembly (/path/to/<multirefdatabase>/refs/*.fa) to add to the yaml')

    _group.add_argument(
        "--species-reference",
        type=as_argparse_type(validate_fasta_file),
        metavar="FASTA",
        default=None,
        required=True,
        help="""(absolute) path to the species reference [required]"""
    )

    _group.add_argument(
        "--clade-reference",
        type=as_argparse_type(validate_fasta_file),
        metavar="FASTA",
        default=None,
        help="""(absolute) path to the clade reference [optional] """
    )

    _group.add_argument(
        '--custom',
        action="store_true",
        help=""" provided reference (and clade) fasta sequence are EXTERIOR, and not included in /path/to/<multirefdatabase> """
    )

    return parser


# ------------------------------------------------------------------------------------- #
# helper functions (mostly for main)
# ------------------------------------------------------------------------------------- #

def get_first_rec_as_row(df:pd.DataFrame) -> NamedTuple:
    """ return the first row from a dataframe as a NamedTuple """
    return next(df.itertuples(index=False, name='Row'))

def get_samtools_stats_summary_stats_dataframe(sam_file:Path) -> pd.DataFrame:
    """ get samtools stats summary numbers in a dataframe """
    summary, comments = get_samtools_stats_summary_numbers(sam_file)
    # reads unmapped:   86248 [A+C]     # <-- in summary
    # view -F 8 -f 4:   53574 [C]
    # view -f 12:       32674 [A]
    # view -f 8 -F 4:   57839 [B]
    # view -f 8:        90513 [A+B]
    df = pd.DataFrame.from_dict(summary,orient="index",columns=['value'])
    df.index.name = "property"
    return df

def match_reference_given_unsorted_sam(args:argparse.Namespace) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    """ """
    # read dataframe containing all fasta accessions linked to speciestags & taxIds
    # TODO: should be read from the <REFERENCE_DATA_DIR>/refs directory, which included the mitochondria
    indices = [ Path(f) for f in glob.glob(os.path.join(REFERENCE_DATA_DIR,"WGS/*.fai")) ]
    taxrecords = [ Path(f) for f in glob.glob(os.path.join(REFERENCE_DATA_DIR,"WGS/*.xml")) if REGEX_PATH_TAX_RECORD.match(f) ]
    if not indices or not taxrecords:
        raise IOError("expected to find *.fai and *.xml files, found %s / %s" % (len(indices),len(taxrecords)))
    df = get_accession_index_dataframe(indices,taxrecords)

    # now get counts from the (potentiall huge) sam file
    counts = get_accession_counter_from_sam_file(args.sam_file)

    #print(df.head())
    #print(df.shape)
    #print(counts.most_common(10))

    # TODO: get "total_reads" from samtools stats SN, containing much more info
    df_stats = get_samtools_stats_summary_stats_dataframe(args.sam_file)
    total_reads = get_paired_end_read_count_from_sam_file(args.sam_file)
    total_mapped = sum(counts.values())

    # connect the counts and the species dataframe
    df['observed'] = df['accession'].map(counts)
    df['observed'] = df['observed'].fillna(0).astype(int)
    df.sort_values(by=['observed'], ascending=[False],inplace=True)
    df.reset_index(inplace=True,drop=True)

    # now group on either/both reference and/or taxid
    # TODO: return both DFs, not just print them
    # TODO: implement argument to choose which of the two DFs to generate and return

    def add_ratio_sort_and_reindex(df:pd.DataFrame) -> None:
        """ add read count ratio's (mapped/allreads), prioritize using sort and reset index """
        df['ratio_mapped'] = round(df.observed / float(total_mapped), 5)
        df['ratio_reads'] = round(df.observed / float(total_reads), 5)
        df.sort_values(by='observed', ascending=False, inplace=True)
        df.reset_index(inplace=True,drop=True)

    grouped_by_reference = df.groupby(['reference', 'taxid'], as_index=False)[['observed']].sum()
    add_ratio_sort_and_reindex(grouped_by_reference)
    grouped_by_taxid = df.groupby(['taxid'], as_index=False)[['observed']].sum()
    add_ratio_sort_and_reindex(grouped_by_taxid)
    return ( grouped_by_reference, grouped_by_taxid, df_stats )


def get_accession_index_dataframe(indices=List[Path],taxrecords=List[Path]) -> pd.DataFrame:
    """ Convert pairs of fasta indices and taxonomy labels to a dataframe of the species assignment search space """
    _dfs = []
    for fai_index,tax_xml in zip(sorted(indices),sorted(taxrecords)):
        _pfx,reference,taxid = REGEX_PATH_TAX_RECORD.match(str(tax_xml)).groups()
        df = pd.read_csv(fai_index, sep='\t', header=None)[[0,1]]
        df.columns = ["accession","length"]
        df['reference'] = reference
        df['taxid'] = int(taxid)
        _dfs.append(df)
    df = pd.concat(_dfs,axis=0)
    # Shurely not the best place, but (some) QC on duplicate fasta accessions should occur
    # In theory >1 fasta can contain the same "some_contig" accession name
    if len(df.accession.unique()) < df.shape[0]:
        raise ValueError("df.accession.unique() isn't unique")
    return df

def get_samtools_stats_summary_numbers(sam_file:Path) -> Tuple[dict,dict]:
    """ get SN summary numbers as dictionary from samtools stats output

    """
    cmd = f""" samtools stats {sam_file} | grep "^SN" | cut -f 2- """
    output = subprocess.check_output(cmd, shell=True, text=True)
    # split in key-value pairs; ignore comments like '# excluding supplementary and secondary reads'
    summary = dict([ tuple(line.split("\t")[0:2]) for line in output.strip().split('\n') ])
    # parse comments, being the 3th column (available for some but not all keys)
    comments = dict([ tuple(line.split("\t")[::2]) for line in output.strip().split('\n') if len(line.split("\t")) == 3 ])
    return (summary,comments)

    #try:
    #    return int(output)
    #except TypeError:
    #    return TypeError(f"expected int, got '{output}'")

def get_paired_end_read_count_from_sam_file(sam_file:Path) -> int:
    """ get the paired end read count from a sam file (thus assuming paired-end Illumina data) """
    cmd = f""" samtools view -c -f 1 {sam_file} """
    output = subprocess.check_output(cmd, shell=True, text=True)
    try:
        return int(output)
    except TypeError:
        return TypeError(f"expected int, got '{output}'")


def get_accession_counter_from_sam_file(sam_file:Path) -> Counter:
    """ Count the number of matches to each accession in a SAM file

    Implementation hard-coded assuming PE-Illumina data!

    # https://broadinstitute.github.io/picard/explain-flags.html

    -F 2048     exclude supplementary alignments
    -F 256      exclude secondary alignments
    -f 2        only reads mapped in proper pairs
    """

    cmd = f""" samtools view -F 2048 -F 256 -f 2 {sam_file} | cut -f 3 """
    p = subprocess.Popen(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, bufsize=1024*1024)

    CHUNK_SIZE = 64 * 1024
    remainder = ""
    counts = Counter()
    with p.stdout:
        while True:
            chunk = p.stdout.read(CHUNK_SIZE)
            if not chunk:
                break
            combined = remainder + chunk
            lines = combined.splitlines(keepends=True)
            if not lines[-1].endswith('\n'):
                remainder = lines.pop()
            else:
                remainder = ""
            # update Counter
            counts.update(line.strip() for line in lines if line.strip())
    # !important! last chunk
    if remainder.strip():
        counts.update([remainder.strip()])
    p.wait()
    return counts


# ------------------------------------------------------------------------------------- #
# output generation functions
# ------------------------------------------------------------------------------------- #

def _generate_yaml_reference_snippet(rec:NamedTuple) -> List[str]:
    """ """
    return [
        # !important! if rec.clade is empty / nan, exclude it in the name
        (f"  name: {rec.species} {rec.clade}").rstrip().rstrip(" nan"),
        f"  accession: {rec.reference}",
        f"  fasta: {rec.fasta}"
        ]

def _write_yaml_from_list(yamltxt:List[str],filename:str,verbose:bool=False) -> None:
    """ write to (systematically layouted) yaml file """
    yamltxt.append("")
    fh = open(filename,'w')
    fh.write("\n".join(yamltxt))
    fh.close()
    if verbose:
        print("\n".join(yamltxt))

# ------------------------------------------------------------------------------------- #
# main match-ref.py function(s)
# ------------------------------------------------------------------------------------- #

def matchref(verbose:bool=True) -> None:
    """ main entry point of this script """
    args = get_parser().parse_args()
    grouped_by_reference, grouped_by_taxid, df_stats = match_reference_given_unsorted_sam(args)

    if verbose:
        print(tabulate(grouped_by_reference.head(10), headers='keys', tablefmt='psql'))
        print(tabulate(grouped_by_taxid.head(10), headers='keys', tablefmt='psql'))
        print(tabulate(df_stats.head(10), headers='keys', tablefmt='psql'))

    # write the dataframe files
    grouped_by_reference.to_csv(args.out_prefix+"-match-ref-reference.tsv",sep='\t')
    grouped_by_taxid.to_csv(args.out_prefix+"-match-ref-taxid.tsv",sep='\t')
    grouped_by_taxid.to_csv(args.out_prefix+"-samtools-stats.tsv",sep='\t')

    # Read dataframe that links genomic/species accessions to reference fasta files
    # !important! at the species level, those defined as cladegroups should result in
    #             species fallback to the "is_primary" accession
    #             The dataframe/tsv (from xlsx) is validated by referencedata.ProvidedSchema
    #             that no inconsistencies are introduced in the dataframe.
    refs_strains = read_reference_species_and_assembly_df()
    refs_species = refs_strains[((refs_strains.is_primary == True) | (refs_strains.is_primary.isnull()))].copy()

    # Get (majority-voted) best matches (species, optionally strain) as NamedTuples
    # From these, link it to the 'refs' dataframe record containing corresponding fasta file names
    # Write this info to <prefix>-references.yml
    # Since the output yml is drop-dead simple, decided to manually write it (not using pyyaml)
    best_matched_species = get_first_rec_as_row(grouped_by_taxid)
    best_species = get_first_rec_as_row(refs_species[refs_species.taxid==best_matched_species.taxid])
    #print(best_matched_species)
    #print(best_species)
    yamltxt = [
        "species:",
        f"  input: {args.sam_file}" ]
    yamltxt.extend( _generate_yaml_reference_snippet(best_species) )

    # !important! Avoid TypeError: boolean value of NA is ambiguous
    if str(best_species.is_primary) == "True":
        best_matched_strain = get_first_rec_as_row(grouped_by_reference)
        try:
            best_strain = get_first_rec_as_row(refs_strains[refs_strains.reference == best_matched_strain.reference])
        except StopIteration:
            # scenario "Alternative Clade I" in reference sheet.
            best_strain = get_first_rec_as_row(refs_strains[refs_strains.taxid == best_matched_strain.taxid])
        #print(best_matched_strain)
        #print(best_strain)
        yamltxt.extend([
            "strain:",
            f"  input: {args.sam_file}",
            ])
        yamltxt.extend(_generate_yaml_reference_snippet(best_strain))

    # write to yml file
    outfname = args.out_prefix + SUFFIX_YAML_MATCHED
    _write_yaml_from_list( yamltxt, outfname )
    sys.stdout.write("EOF: written %s\n" % outfname)


def generate_forced_reference_file() -> None:
    """ generate a "forced" reference/strain yaml file which overrules match-ref.py vanilla *-references.yml file """
    parser = get_forced_parser()
    args = parser.parse_args()

    if args.custom:
        # exterior provide fasta sequences.
        # - Mimic Dummy "record" for attribute mapping
        # - add exterior: true
        class Dummy:
            reference = "NA"
            species = "NA"
            clade = "" #!important!
        Dummy.fasta = os.path.basename(args.species_reference)
        yamltxt = [ "species:" ]
        yamltxt.extend( _generate_yaml_reference_snippet(Dummy) )
        yamltxt.extend(["  exterior: true"])
        if args.clade_reference:
            Dummy.fasta = os.path.basename(args.clade_reference)
            yamltxt.extend(["clade:"])
            yamltxt.extend(_generate_yaml_reference_snippet(Dummy))
            yamltxt.extend(["  exterior: true"])

    else:
        # fasta sequences should be in the multireference database. Check this
        msg_not_in_refdb = ( "fasta not occurring in %s " % REFERENCE_DATA_DIR ) + "[%s]"
        refs_strains = read_reference_species_and_assembly_df()
        refs_species = refs_strains[((refs_strains.is_primary == True) | (refs_strains.is_primary.isnull()))].copy()

        try:
            species = get_first_rec_as_row(refs_species[refs_species.fasta == os.path.basename(args.species_reference)])
        except StopIteration:
            parser.error(msg_not_in_refdb % args.species_reference)
        yamltxt = [ "species:" ]
        yamltxt.extend( _generate_yaml_reference_snippet(species) )
        if args.clade_reference:
            try:
                clade = get_first_rec_as_row(refs_strains[refs_strains.fasta == os.path.basename(args.clade_reference)])
            except StopIteration:
                parser.error(msg_not_in_refdb % args.clade_reference)
            yamltxt.extend(["clade:"])
            yamltxt.extend(_generate_yaml_reference_snippet(clade))

    # write to yml file
    outfname = args.out_prefix + SUFFIX_YAML_FORCED
    _write_yaml_from_list( yamltxt, outfname)
    sys.stdout.write("EOF: written %s\n" % outfname)

def main():
    if "--species-reference" in sys.argv:
        # switch to generating a forced/bypass yaml file
        generate_forced_reference_file()
    else:
        # Your original main logic goes here
        matchref()

if __name__ == "__main__":
    main()
