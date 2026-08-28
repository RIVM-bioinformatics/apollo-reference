#!/usr/bin/env python3

import os
import sys
import hashlib
from typing import Optional
import pandas as pd
import pandera.pandas as pa
from pandera.typing import Series
from pathlib import Path
from typing import NewType, Dict, Any, Union
from rich.console import Console
from rich.markdown import Markdown

# Package Imports
try:
    from configuration import ASSEMBLY_REFERENCE_TSV_BASE, SPECIES_REFERENCE_TSV_BASE
except ModuleNotFoundError:
    from .configuration import ASSEMBLY_REFERENCE_TSV_BASE, SPECIES_REFERENCE_TSV_BASE

# Customized TypeHints
VerifiedFile = NewType("VerifiedFile", Path)

## Define
#refdata_schema = pa.DataFrameSchema({
#    "speciesTaxId": pa.Column(int, checks=pa.Check.greater_than(0)),
#    "Email": pa.Column(str, checks=pa.Check.str_matches(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$")),
#    "Status": pa.Column(str, checks=pa.Check.isin(["Actief", "Inactief"])),
#})

accession_regex_ncbi=r"GC[AF]_\d{9}\.\d"
accession_regex_rivm=r"RIVM\d{6}"
accession_regex = r"^("+"|".join([accession_regex_ncbi,accession_regex_rivm])+")$"
accession_regex_MT = accession_regex[0:-2]+"|included|)$"

class ProvidedSchema(pa.DataFrameModel):
    """ species configuration dataframe (external), "reference-tsv", as provided by the user """
    reference: Series[str] = pa.Field(regex=accession_regex)
    mitochondrion: Series[str] = pa.Field(regex=accession_regex_MT,nullable=True)
    datasource: Series[str] = pa.Field(isin=["NCBI", "RIVM"])
    ploidy: Series[str] = pa.Field(isin=["haploid", "diploid","haploid/diploid"],nullable=True)
    cladegroup: Series[str] = pa.Field(nullable=True)
    is_primary: Optional[Series["Int64"]] = pa.Field(isin=[0,1], nullable=True, coerce=True)
    ignore: Optional[Series["Int64"]] = pa.Field(isin=[0,1], nullable=True, coerce=True)
    # all the biology-friendly fields
    category: Series[str] = pa.Field(str_length={"max_value": 200})
    literature: Series[str] = pa.Field(str_length={"max_value": 200},nullable=True)
    species: Series[str] = pa.Field(str_length={"max_value": 150})
    clade: Series[str] = pa.Field(str_length={"max_value": 50},nullable=True)
    strain: Series[str] = pa.Field(str_length={"max_value": 50},nullable=True) # preferably require this?
    recent_name: Series[str] = pa.Field(str_length={"max_value": 150},nullable=True)

    class Config:
        coerce = True  # Belangrijk: pandas zet int-kolommen met nulls vaak om naar floats

    @pa.dataframe_check
    def check_both_none_or_neither(cls, df):
        """ is_primary and cladegroup can be None only if both are None """
        #both_none = df["is_primary"].isna() & df["cladegroup"].isna()
        #neither_none = df["is_primary"].notna() & df["cladegroup"].notna()
        #return both_none | neither_none
        return (df["is_primary"].isna() == df["cladegroup"].isna()).all()

    @pa.check("is_primary", groupby="cladegroup", name="cladegroup_master_informants_requirement")
    def check_per_category_logic(cls, is_primary_grouped: Dict[Any, pd.Series]) -> bool:
        """
        When groupby is provided, Pandera passes a dict: {group_name: sub_series}.
        Each category must have exactly one '1' and at least one '0'.
        """
        for cladegroup, group in is_primary_grouped.items():
            # Skip rows where cladegroup was None (Pandera handles this, but good to be safe)
            if cladegroup is None or pd.isna(cladegroup):
                continue
            count_one = (group == 1).sum()
            count_zero = (group == 0).sum()
            if not (count_one == 1 and count_zero >= 1):
                return False
        return True

    #@pa.check("column3")
    #def custom_check(cls, series: pd.Series) -> pd.Series:
    #    return series.str.len() == 1

class DerivedSchema(pa.DataFrameModel):
    """ species configuration dataframe "reference_assembly_data.tsv" (internal), containing data on the actual assemblies """
    reference: Series[str] = pa.Field(regex=accession_regex)
    taxid: int = pa.Field(ge=0)
    MT_num: int = pa.Field(ge=0,nullable=True)
    MT_assembly: str = pa.Field(str_length={"max_value": 200},nullable=True)
    MT_accession: str = pa.Field(str_length={"max_value": 50},nullable=True) # maybe comma-sep in case fragmented?
    MT_length: int = pa.Field(ge=0,nullable=True)
    MT_source: str = pa.Field(isin=["included", "external", "provided"],nullable=True)
    fasta: str = pa.Field(str_length={"max_value": 200})

def validate_df_using_pa(df:pd.DataFrame, schema:pa.DataFrameSchema, lazy:bool=True, dryrun:bool=False) -> bool:
    """ """
    try:
        schema.validate(df, lazy=lazy) # lazy=True reports all errors, not just the first
        return True
    except pa.errors.SchemaErrors as err:
        if dryrun:
            print(f"Validation failed:\n{err.failure_cases}")
            return False
        else:
            raise pa.errors.SchemaErrors(err.failure_cases)


def get_species_index_hash_from_df(df:pd.DataFrame) -> str:
    """
    Generate unique reference accession dataset hash, that describes the exact content of the database

    !important! mind that this should (not DRY ...) exactly mimic the function
                generate_reference_accession_dataset_UID()
                in workflow/scripts/build-helpers.sh
    """
    validate_df_using_pa(df,DerivedSchema)
    accessions = list(df[['reference','MT_accession']].sort_values(['reference']).itertuples(index=False, name=None))
    accessions = ",".join([ "%s,%s" % t for t in accessions ]).replace("nan,","")
    return hashlib.md5(accessions.encode()).hexdigest()

def validate_reference_dataset(dbpath:Union[Path|str],df:pd.DataFrame=None,default_df_name:str=ASSEMBLY_REFERENCE_TSV_BASE) -> Path:
    """ validate that the dataframe matches the provided reference database dir, and that this dir has the proper blueprint """
    dbpath = str(dbpath)
    # Check all required subdirectories; SRR is an optional sub-directory!
    for subdir in ('','WGS','refs','mmidx','NUCCORE',):
        dirpath = os.path.join(dbpath,subdir)
        if not os.path.isdir(dirpath):
            raise FileNotFoundError(dirpath)
        elif subdir != '' and len(os.listdir(dirpath)) == 0:
            raise FileNotFoundError("no files found in %s" % dirpath)

    if type(df) == type(None) and os.path.isfile(os.path.join(dbpath,default_df_name)):
        # Dataframe is expected (vanilla) to be IN the directory;
        # !important! delayed import prevents circular import
        from .dataframes import read_reference_assembly_df
        dfpath = os.path.join(dbpath,default_df_name)
        df = read_reference_assembly_df(dfpath)
    elif type(df) == type(None):
        # Fallback to (nearly requiredly present...) dataframe impossible since file is not there
        msg = "No assembly dataframe provided; thus can't validate exact files being present"
        raise ValueError(msg)

    # check if a correct dataframe format was provided
    validate_df_using_pa(df,DerivedSchema)

    # validate that all 'fasta' constructed references are at their designated location (in /refs)
    for fasta in df.fasta.tolist():
        absfasta = os.path.join(dbpath,'refs',fasta)
        if not os.path.isfile(absfasta):
            raise FileNotFoundError(absfasta)

    # validate that the identify-species mmidx corresponds to an already existing one
    mmidx = get_identify_species_mmidx_relpath(df)
    mmidx = os.path.join(dbpath, 'mmidx', mmidx)
    if not os.path.isfile(mmidx):
        raise FileNotFoundError(mmidx)
    return Path(dbpath)


def validate_reference_dataset_multiclade_requirements(dbpath:Union[Path|str],df:pd.DataFrame=None,default_df_name:str=SPECIES_REFERENCE_TSV_BASE) -> bool:
    """ validate that files needed in the pipeline for multiclade-analyses indicated in the dataframe, exist """
    dbpath = str(dbpath)

    fullpath_default_df_name = os.path.join(dbpath,default_df_name)
    if type(df) == type(None) and os.path.isfile(fullpath_default_df_name):
        # Dataframe is expected (vanilla) to be IN the directory;
        # !important! delayed import prevents circular import
        from .dataframes import read_reference_species_df
        df = read_reference_species_df(fullpath_default_df_name)
    elif type(df) == type(None):
        # Fallback to (nearly requiredly present...) dataframe impossible since file is not there
        print(fullpath_default_df_name)
        msg = "No assembly dataframe provided; thus can't validate is multiclade masking files being present"
        raise ValueError(msg)

    # check if a correct dataframe format was provided
    # TODO: change DerivedSchema and the TSV to include 'is_primary' and 'cladegroup'
    #validate_df_using_pa(df,DerivedSchema)
    validate_df_using_pa(df,ProvidedSchema)

    # TODO: need to get stated somewhere in config, not hardcoded. Now:
    # - apollo-mapping having to know these exact filename too
    # - bash script that generates this file has to know filename too
    # - and files currently (and temporarily!!) reside in apollo-mapping, not *-reference ...
    # - ... since their generation is not fully automated yet.
    blacklist_filenames = [
        os.path.join(dbpath,"multiclade","cauris-GCA_002759435.3-vs-WGA-blacklist.bed"),
        os.path.join(dbpath, "multiclade", "cauris-GCA_002759435.3-vs-fastq-blacklist.bed"),
        os.path.join(dbpath, "multiclade", "cauris-GCA_002759435.3-blacklist.bed"),
        ]
    for fname in blacklist_filenames:
        if not os.path.isfile(fname):
            # Required blacklist file not there.
            # - Future solution is that this is "impossible" since apollo-reference made its own blacklist(s) automatically
            # - Current situation is that (semi-manually) generated files for C.auris are stored in apollo-mapping/data
            # - So, non-present file most likely means a newly downloaded and/or externally downloaded multireference-dataset ...
            # - for which the file's aren't copied yet (see apollo-mapping/data/README.md)
            if os.path.basename(os.getcwd()) == "apollo-mapping":
                MARKDOWN_FILE = os.path.join(os.getcwd(),"data","README.md")
                error_msg = f"""## This is the exception that is about to get raised:\nraise FileNotFoundError({fname})\n\n"""
                md = Markdown(error_msg + open(MARKDOWN_FILE).read())
                console = Console()
                console.print(md)
                sys.exit()

            # vanilla / any other situation. Required file is missing
            raise FileNotFoundError(fname)

    # Okay(ish) for now. In the future it will be better to fully control
    # the existance of all of the components that make up the blacklist files.from
    return True

    # --------------------------------------------------------------------------------------------------- #
    # Since there is at this moment:
    # - only 1 multiclade species supported (C.auris)
    # - and the "multiclade" concept is decided to be an expert feature
    # - and the generation of these blacklist files is not 100% automated,
    # - since it was decided thyis was "out-of-scope" of this batch of work on apollo-mapping,
    # For now these blacklists are stated in the repository itself.
    # Users have to, when wanting to use this, copy the files to the designated location
    # (in their local multispecies-database), and from their on have a fully functional database.
    # --------------------------------------------------------------------------------------------------- #
    # For how the eventual validation could look like, below code snippets are fairly representative
    # Check if all the required files for a succesful multiclade-SNP-analyses are present
    # TODO: code below needs to get refactored to the following logics:
    #       - there are pairwise comparison files of SRR vs central reference (assembly)
    #       - there are pairwise comparison files of assembly vs central reference (assembly)
    #       - there are "merged" files per accession/assembly/SRR combination
    #       momentarily it only checks a subset, and only sticks to strict "GCA" accession file naming

    # Suffixes of files required to be present.
    # Example: multiclade/GCA_003013715.2__vs__GCA_002759435.3-sofclippedregions.bed
    # TODO: DRY define these in a yaml config;
    #       This allows the filename to be written to correspond to the filename to be validated here
    # TODO: Eventually it's best practice to truely validate all these files (being present) individually;
    #       For now, this will be simplified to only checking if the overall BED blacklist is there
    required_suffixes = [
        'sofclippedregions.bed',
        'segmental-deletions.bed',      # TODO: of at least size Xnt; where to decide and/or filter on this?
        'segmental-insertions.bed',     # TODO: of at least size Xnt; where to decide and/or filter on this?
        'segmental-duplications.bed',   # TODO: of at least size Xnt; where to decide and/or filter on this?
        #'hypervariable.bed'
        'noncovered.bed'
    ]

    fdf = df[df.is_primary.notnull()]
    for idx,row in fdf[fdf.is_primary == 0].iterrows():
        _, primary = next(fdf[((fdf.cladegroup==row['cladegroup']) & (fdf.is_primary==1))].iterrows())
        central_accession = primary['reference']
        for suffix in required_suffixes:
            clade_accession = row['reference']
            fname = os.path.join(dbpath, 'multiclade', f"{clade_accession}__vs__{central_accession}-{suffix}")
            if not os.path.isfile(fname):
                raise FileNotFoundError(fname)
    return True

def get_identify_species_mmidx_relpath(df:pd.DataFrame=None) -> str:
    """ generate the apollo-species-refs.<md5hash>.mmidx from the provided dataframe

    Dataframe should be the one derived from:
        - read_reference_assembly_df
        - read_reference_species_and_assembly_df
    """
    # check if the correct dataframe format was provided
    validate_df_using_pa(df, DerivedSchema)
    md5hash = get_species_index_hash_from_df(df)
    return 'apollo-species-refs.%s.mmidx' % md5hash

"""
irods_collection_name_format_<hash>
    genomes
        + raw, putatively including mitochondria
    mitochondria
        + raw, explicitly downloaded mitochondria
    references
        + reference genome including corresponding mitochondria
    mmidx
        1 identify_species_minimap2_index_<hash>.mmidx
    SRR
        ~ downloaded fastq corresponding to reference genomes (optional, for testing)
    sam
        ~ mapped (SRR) reads to mmidx "match-ref" (optional, for testing)
"""

# TODO: Eventual rewrite location of reference data to IRODS collection
#       This means workflow should be rewritten to be compatible to either
#        - mount or share folder
#        - IRODS collection
#       This is being (naively ...) prepared in ReferenceDataHandler class
#       Issue here is that, once (e.g.) apollo_mapping.py kicks in:
#         - it can't validate "referencedata" in case it's not on a mount
#         - can't even provide --help, since it needs TSV files linked to referencedata

class ReferenceDataHandler:
    """ ReferenceData handler that is basically "read-only", and can only get access to existing ReferenceData builds  """
    def __init__(self, df:pd.DataFrame):
        self.df = df

    def validate(self):
        return validate_df_using_pa(self.df,Schema)

    def generate_hash(self) -> str:
        """ generate a persistent hash from all the accessions in the dataset """
        data = "=".join(self.df.accession.to_list(sorted=True))
        return hashlib.sha256(data.encode('utf-8')).hexdigest()

    def get_local_path_to_xxxx(self) -> VerifiedFile:
        """ """
        raise NotImplementedError

    def get_minimap2_index_name(self) -> str:
        """ """
        return "identify_species_minimap2_index_%s.mmidx" % self.generate_hash()

    def require_minimap2_index(self)-> VerifiedFile:
        """ Builds and/or retrieves the required minimap2 reference species index

        returns the full path to the crucially required existing & accessible minimap2 index
        """
        raise NotImplementedError

    def list_reference_genomes(self,include_mitochondria:bool=False) -> pd.DataFrame:
        """ return a dataframe containing accession codes and full paths to reference genomes """
        raise NotImplementedError

class ReferenceDataBuilder(ReferenceDataHandler):
    """ ReferenceData handler that can as well build new indices """

    def validate(self):
        super(self).validate()
        # TODO: validate that downloaded genomes are fasta?
        # TODO: validate that .....

    def build_minimap2_index(self) -> VerifiedFile:
        """ generates a new minimap2 index containing all references """
        indexname = "identify_species_minimap2_index_%s.mmidx" % self.generate_hash()
        raise NotImplementedError

    def download_accession(self,accession) -> bool:
        """ downloads a new stated reference WGS from NCBI using datasets """
        raise NotImplementedError

class ReferenceDataPublisher(ReferenceDataHandler):
    """ ReferenceData handler that can publish a new build to a place where it can be retrieved from at a later point in time """
    def validate(self):
        super(self).validate()
        # TODO: validate that all data to publish is on the file system
        # TODO: validate that location where to publish to is accessible

    def publish_minimap2_index(self) -> VerifiedFile:
        """ make sure to "publish" the required minimap2 index on your filesystem """
        # download from iRODs
        # build de-novo
        raise NotImplementedError
