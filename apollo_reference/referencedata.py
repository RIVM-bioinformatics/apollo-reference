#!/usr/bin/env python3

import os
import hashlib
from subprocess import Popen, PIPE
from typing import Optional
import pandas as pd
# import pandera.pandas as pa was syntax in older versions!
import pandera as pa
from pandera.typing import Series
from pathlib import Path
from typing import NewType, Dict, Any, Union

# Customized TypeHints
VerifiedFile = NewType("VerifiedFile", Path)

## Define
#refdata_schema = pa.DataFrameSchema({
#    "speciesTaxId": pa.Column(int, checks=pa.Check.greater_than(0)),
#    "Email": pa.Column(str, checks=pa.Check.str_matches(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$")),
#    "Status": pa.Column(str, checks=pa.Check.isin(["Actief", "Inactief"])),
#})

accession_regex_ncbi="GC[AF]_\d{9}\.\d"
accession_regex_rivm="RIVM\d{6}"
accession_regex = r"^("+"|".join([accession_regex_ncbi,accession_regex_rivm])+")$"
accession_regex_MT = accession_regex[0:-2]+"|included|)$"

class ProvidedSchema(pa.DataFrameModel):
    """ species configuration dataframe (external), "reference-tsv", as provided by the user """
    reference: Series[str] = pa.Field(regex=accession_regex)
    mitochondrion: Series[str] = pa.Field(regex=accession_regex_MT,nullable=True)
    datasource: Series[str] = pa.Field(isin=["NCBI", "RIVM"])
    ploidy: Series[str] = pa.Field(isin=["haploid", "diploid"],nullable=True)
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

def validate_reference_dataset(dbpath:Union[Path|str],df:pd.DataFrame=None) -> bool:
    """ validate that the provided reference database dir has the blueprint of it """
    dbpath = str(dbpath)
    # SRR is an optional sub-directory!
    for subdir in ('','WGS','refs','mmidx','NUCCORE',):
        dirpath = os.path.join(dbpath,subdir)
        if not os.path.isdir(dirpath):
            raise FileNotFoundError(dirpath)
        elif subdir != '' and len(os.listdir(dirpath)) == 0:
            raise FileNotFoundError("no files found in %s" % dirpath)

    if type(df) == type(None):
        # no assembly dataframe provided; thus can't check exact files being present
        return True

    # check if the correct dataframe format was provided
    validate_df_using_pa(df,DerivedSchema)

    # validate that all 'fasta' constructed references are at their designated spot
    for fasta in df.fasta.tolist():
        absfasta = os.path.join(dbpath,'refs',fasta)
        if not os.path.isfile(absfasta):
            raise FileNotFoundError(absfasta)

    # validate that the identify-species mmidx corresponds to an already existing one
    mmidx = get_identify_species_mmidx_relpath(df)
    mmidx = os.path.join(dbpath, 'mmidx', mmidx)
    if not os.path.isfile(mmidx):
        raise FileNotFoundError(mmidx)
    return True

def get_identify_species_mmidx_relpath(df:pd.DataFrame=None) -> str:
    """ generate the apollo-species-refs.<md5hash>.mmidx from the provided dataframe """
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
