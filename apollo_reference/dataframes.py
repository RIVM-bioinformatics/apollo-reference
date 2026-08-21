from __future__ import annotations

# Python Imports
#import os
import re
import pandas as pd
from typing import TYPE_CHECKING

# Package Imports
try:
    from referencedata import ProvidedSchema
    from configuration import SPECIES_REFERENCE_TSV, ASSEMBLY_REFERENCE_TSV, FASTQ_REFERENCE_TSV
except ModuleNotFoundError:
    from .referencedata import ProvidedSchema
    from .configuration import SPECIES_REFERENCE_TSV, ASSEMBLY_REFERENCE_TSV, FASTQ_REFERENCE_TSV

if TYPE_CHECKING:
    # doesn't happen during runtime
    import pandera.pandas as pa

def read_reference_species_df(fname:str=SPECIES_REFERENCE_TSV,schema:pa.DataFrameModel=ProvidedSchema) -> pd.DataFrame:
    """ read the/a reference species table and convert human-readable to machine-readable conventions """
    df = pd.read_csv(fname, sep="\t")
    df._metadata.append("tsv_source")
    df.tsv_source = fname
    # strip all exterior whitespace from text
    # !important! panda version >= 2.1.0 renamed to map
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    df = df[df["ignore"].ne(1)]
    df['is_primary'] = df['is_primary'].map({1.0: True, 0.0: False}).astype('boolean')
    pa_column_names = list(schema.to_schema().columns.keys())
    df_column_names = df.columns.tolist()
    renames = {}
    for colname in df_column_names:
        if colname not in pa_column_names and colname.lower() in pa_column_names:
            # example: Datasource -> datasource
            renames[colname] = colname.lower()
        elif colname not in pa_column_names and colname.split()[0].lower() in pa_column_names:
            # example: Species name -> species
            renames[colname] = colname.split()[0].lower()
        elif colname not in pa_column_names and re.search("^(recent|current) name",colname.lower()):
            # example: Current name (Kidd et al. 2023) -> recent_name
            renames[colname] = "recent_name"
    if renames:
        df.rename(columns=renames, inplace=True)
    return df

def read_reference_assembly_df(fname:str=ASSEMBLY_REFERENCE_TSV) -> pd.DataFrame:
    """ read the machine-generated reference assembly table (thus unvalidated) """
    df = pd.read_csv(fname, sep="\t")
    df._metadata.append("tsv_source")
    df.tsv_source = fname
    df['MT_num'] = df['MT_num'].astype('Int64')
    df['MT_length'] = df['MT_length'].astype('Int64')
    return df

def read_reference_species_and_assembly_df() -> pd.DataFrame:
    """ TODO: read from input files as arguments """
    df = read_reference_species_df()
    ProvidedSchema.validate(df, lazy=True)
    df.drop(columns=['ignore', 'literature', 'recent_name'], inplace=True)
    df.drop(columns=['datasource', 'category', 'ploidy'], inplace=True)
    # pitch in MT assembly information
    dfasm = read_reference_assembly_df()
    df = pd.merge(df, dfasm, on='reference', how='left')
    df.tsv_source = (df.tsv_source,dfasm.tsv_source)
    df['taxid'] = df['taxid'].astype('Int64')
    df = df[~df.taxid.isnull()]
    return df
