#!/usr/bin/env python3
""" summarize all the downloaded assembly data in one single tab-delimited dataframe """

# Python Imports
import pandas as pd
import os
import glob
import sys
from tabulate import tabulate
from pathlib import Path

def generate_dataframe(REFDATADIR:Path) -> pd.DataFrame:
    """ Convert all downloaded assembly data, genomic and mitochondria, into a single, summarizing pandas dataframe

    Example how this dataframe looks like (in csvlook -t -I on the command line):

    | reference       | taxid  | MT_num | MT_assembly                                     | MT_accession   | MT_length | MT_source | fasta                                        |
    | --------------- | ------ | ------ | ----------------------------------------------- | -------------- | --------- | --------- | -------------------------------------------- |
    | GCA_003054445.1 | 4909   | 1      | GCA_003054445.1_ASM305444v1_genomic.fna         | CP028778.1     | 51340     | included  | GCA_003054445.1__included__CP028778.1.fa     |
    | GCA_010111755.1 | 5478   | 15     | GCA_010111755.1.mito.3062092141.fa              | CP130941.1     | 20062     | external  | GCA_010111755.1__external__CP130941.1.fa     |
    | GCA_004931855.1 | 52247  | 1      | GCA_004931855.1.mito.3111511661.fa              | BK073862.1     | 40077     | external  | GCA_004931855.1__external__BK073862.1.fa     |
    | GCA_030570855.1 | 223819 |        |                                                 |                |           |           | GCA_030570855.1__nd__nd.fa                   |
    | GCA_000182965.3 | 237561 | 2      | GCA_000182965.3.mito.266630551.fa               | NC_002653.1    | 40420     | external  | GCA_000182965.3__external__NC_002653.1.fa    |
    | GCA_030569875.1 | 273131 | 3      | GCA_030569875.1.mito.1314948788.fa              | NC_036380.1    | 46138     | external  | GCA_030569875.1__external__NC_036380.1.fa    |
    | GCA_000006445.2 | 284592 |        |                                                 |                |           |           | GCA_000006445.2__nd__nd.fa                   |

    """
    sys.stderr.write("# Generating reference assembly dataframe: "+REFDATADIR+"\n")
    sys.stderr.write("# Generating reference assembly dataframe: "+REFDATADIR+"\n")
    sys.stderr.write("# Generating reference assembly dataframe: "+REFDATADIR+"\n")
    # read files that link reference accession to taxId and cast into dataframe
    data = []
    for xml in glob.glob(os.path.join(REFDATADIR, "WGS", "*.*.xml")):
        _parts = Path(xml).stem.split(".")
        accession = ".".join(_parts[0:-1])
        taxId = _parts[-1]
        data.append((accession, int(taxId)))
    colnames = ['reference', 'taxid']
    df = pd.DataFrame(data, columns=colnames)
    df.sort_values(list(reversed(colnames)), inplace=True)

    def read_mt_df(localfname:str) -> pd.DataFrame:
        """ DRY helper function to read provided.*, included.* and selected.mitochondria.tsv files """
        colnames_mt = "reference num assembly accession length source".split()
        df_MT = pd.read_csv(os.path.join(REFDATADIR, localfname), header=None, sep="\t")
        df_MT['source'] = localfname.split(".")[0]
        df_MT.columns = colnames_mt
        return df_MT

    df_MT_inc = read_mt_df("included.mitochondria.tsv")
    df_MT_prv = read_mt_df("provided.mitochondria.tsv")
    df_MT_ext = read_mt_df("selected.mitochondria.tsv")
    df_MT_ext.source = 'external'

    # this is what we want to achieve:
    # prioritize to those with provided > included > external MT in case multiple available
    # realize things can bumb here nastily in case provided + internal clash
    # scenario 1: provided & internal, but the same MT accession
    df_MT = pd.concat([df_MT_prv,df_MT_inc], axis=0)
    df_MT.sort_values(['reference','source'], ascending=[True,False],inplace=True)
    df_MT_provided_is_included = df_MT[df_MT.duplicated(subset=['reference','accession'], keep=False)]
    if not df_MT_provided_is_included.empty:
        # remove from df_MT_prv (since MT accession is already included!)
        idx_to_drop = df_MT_prv[df_MT_prv['reference'].isin(df_MT_provided_is_included.reference.unique())].index
        df_MT_prv.drop(index=idx_to_drop, inplace=True)

    # scenario 2: provided & internal, but DISTINCT MT accession
    df_MT = pd.concat([df_MT_prv,df_MT_inc], axis=0)
    uniq_counts = df_MT.groupby('reference')['accession'].transform('nunique')
    df_MT_clash = df_MT[uniq_counts>1]
    if not df_MT_clash.empty:
        sys.stderr.write("# WARNING: provided AND included MT clash:\n" + str(df_MT_clash.head(40)) + "\n")
        # remove from df_MT_prv (better to rely on what's provided in NCBI submissions ...)
        idx_to_drop = df_MT_prv[df_MT_prv['reference'].isin(df_MT_clash.reference.unique())].index
        df_MT_prv.drop(index=idx_to_drop, inplace=True)

    # finally ready now to:
    # prioritize to those with provided > included > external MT in case multiple available
    df_MT = pd.concat([df_MT_prv,df_MT_inc, df_MT_ext], axis=0)
    df_MT.sort_values(['reference','source'], ascending=[True,False],inplace=True)
    df_MT = df_MT.drop_duplicates(subset=['reference'], keep='first')

    # merge MT dataframe into genomic dataframe and to final type patches
    df_MT.columns = [ "MT_"+col for col in df_MT.columns.tolist() ]
    df = pd.merge(df, df_MT, left_on='reference', right_on='MT_reference', how='left')
    df.drop(['MT_reference'], axis=1, inplace=True)
    df['MT_num'] = df['MT_num'].astype('Int64')
    df['MT_length'] = df['MT_length'].astype('Int64')
    # eventual actual combined reference name
    df['fasta'] = df['reference'].str.cat([df['MT_source'].fillna('nd'), df['MT_accession'].fillna('nd')], sep='__') + ".fa"
    df.sort_values(list(reversed(colnames)), inplace=True)
    return df

if __name__ == "__main__":
    # (global) variables
    try:
        REFDATADIR = sys.argv[1]
        if not os.path.isdir(REFDATADIR):
            raise IOError
    except:
        help_text = "%s [/path/to/multiref/dir]"
        print(help_text % sys.argv[0])
        sys.exit()

    df = generate_dataframe(REFDATADIR)
    if '--tabulate' in sys.argv:
        # not documented, developer flag only
        print(tabulate(df, headers='keys', showindex=False, tablefmt='psql'))
    else:
        print(df.to_csv(index=False,sep="\t").rstrip())