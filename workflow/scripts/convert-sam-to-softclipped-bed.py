#!/usr/bin/env python3

import os
import io
import sys
from pathlib import Path
from typing import Tuple
import pandas as pd
import numpy as np
import subprocess
import argparse

def read_genome_dict_from_sambam_file(bam_file:Path) -> dict:
    """ """
    cmd = ["samtools", "view", "-H", bam_file]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, _ = p.communicate()

    genome_dict = {}
    for line in stdout.splitlines():
        if line.startswith("@SQ"):
            parts = line.split("\t")
            chrom = [p.replace("SN:", "") for p in parts if p.startswith("SN:")][0]
            length = [int(p.replace("LN:", "")) for p in parts if p.startswith("LN:")][0]
            genome_dict[chrom] = length

    return genome_dict

def naive_is_headered_sambam_file(parser:argparse.ArgumentParser,arg:str,extensions=['.bam','.sam']) -> Path:
    """ Naively test if the provided argument is a headered SAM/BAM file """
    cmd = f"samtools view -H {arg}"
    if not os.path.isfile(arg):
        parser.error(f"not a file: {arg}")
    elif not os.path.splitext(arg)[-1] in extensions:
        parser.error(f"not a SAM/BAM file: {arg}")
    elif subprocess.run(cmd, shell=True, capture_output=True, text=True).stderr:
        parser.error(f"not a (headered) SAM/BAM file: {arg}")
    elif not subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout:
        parser.error(f"not a headered SAM/BAM file: {arg}")
    return Path(arg)

def naive_is_headered_sam_file(parser:argparse.ArgumentParser,arg:str) -> Path:
    """ Naively test if the provided argument is a headered SAM file """
    return naive_is_headered_sambam_file(arg,extensions=['.sam'])

def naive_is_headered_bam_file(parser:argparse.ArgumentParser,arg:str) -> Path:
    """ Naively test if the provided argument is a headered BAM file """
    return naive_is_headered_sambam_file(arg,extensions=['.bam'])

def get_parser() -> argparse.ArgumentParser:
    """ """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "sambam_file",
        type=lambda x: naive_is_headered_sambam_file(parser, x),
        metavar="SAM/BAM",
        help="(unsorted) headered SAM/BAM file",
    )
    return parser

def sam2softclippedbed(sambam:Path) -> pd.DataFrame:
    """ """
    genome_size_dict = read_genome_dict_from_sambam_file(sambam)

    cmd = """ samtools view -e 'cigar =~ "S"' -h %s | samtools view -Sbh | bedtools bamtobed -i - -cigar """ % sambam
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = p.communicate()
    del(errs)
    df = pd.read_csv(io.BytesIO(outs), sep='\t', header=None)
    FINAL_COLORDER = ['chrom','sta','end','name','score','strand','cigar']
    df.columns = FINAL_COLORDER
    df['clip5p'] = df['cigar'].str.extract(r'^(\d+)S').astype(float)
    df['clip3p'] = df['cigar'].str.extract(r'(\d+)S$').astype(float)

    def get_clipped_coords(row:pd.Series) -> Tuple[int,int]:
        """ convert clipping lengths to all-clipped coordinate sta,end ranges """
        if (row['strand'],pd.isna(row['clip5p'])) == ('+',True):
            return ( int(row['end']), int(row['end'])+int(row['clip3p']) )
        elif (row['strand'], pd.isna(row['clip5p'])) == ('-', True):
            return ( int(row['sta'])-int(row['clip3p']), int(row['sta']) )
        elif (row['strand'], pd.isna(row['clip3p'])) == ('+', True):
            return ( int(row['sta'])-int(row['clip5p']), int(row['sta']) )
        elif (row['strand'], pd.isna(row['clip3p'])) == ('-', True):
            return ( int(row['end']), int(row['end'])+int(row['clip5p']) )
        else:
            # !important! mapping clipped on both sides,
            # so will end up as two all-clipped coordinate sta,end ranges
            return (0,0)

    # append original coordinates as concatenated string (track-and-trace input to output)
    df['orig'] = df['sta'].astype(str) + ',' + df['end'].astype(str)
    # now deal with dual-side clippings. These will yield original rows to get duplicated!
    ext5p = df[ ~df.clip5p.isna() & ~df.clip3p.isna() ].copy()
    ext3p = ext5p.copy()
    df = df[ df.clip5p.isna() | df.clip3p.isna() ]
    ext5p.clip3p = pd.NA
    ext3p.clip5p = pd.NA
    for _df in (df, ext5p, ext3p):
        _df[['clip_sta', 'clip_end']] = _df.apply(get_clipped_coords, axis=1, result_type='expand')

    if False:
        for _df in (ext5p, ext3p):
            #print(_df.shape)
            #print(_df.head(20))
            #print(_df[_df.clip_sta == _df.clip_end].shape)
            #print(_df[_df.clip_sta == _df.clip_end].head(20))
            #copied = _df.copy()
            _df['clip_sta'] = _df['clip_sta'].clip(lower=0)
            _df['clip_end'] = np.minimum(_df['clip_end'], _df['chrom'].map(genome_size_dict))
            # remove entries with 0,0 or end,end coordinates (clipping at the very exterior of the sequences)
            _df  = _df[_df.clip_sta != _df.clip_end]
            #print(_df[_df.clip_sta == _df.clip_end].shape)
            #print(_df[_df.clip_sta == _df.clip_end].head(20))
            #print(copied[copied.index.isin(_df[_df.clip_sta == _df.clip_end].index.tolist())].head(20))
        #print(ext5p.head(20))
        #print(ext3p.head(20))
        sys.exit()

    # finalize: concat, set 0<=..<=end coordinates, sort and make BED6+cigar+orig
    df = pd.concat([df,ext5p,ext3p])
    df['clip_sta'] = df['clip_sta'].clip(lower=0)
    df['clip_end'] = np.minimum(df['clip_end'], df['chrom'].map(genome_size_dict))
    # remove entries with 0,0 or end,end coordinates (clipping at the very exterior of the sequences)
    df = df[df.clip_sta != df.clip_end]
    df.sort_values(by=['chrom','clip_sta','clip_end'],ascending=[True,True,False], inplace=True)

    # refactor to final dataframe
    df.drop(columns=['sta','end','clip5p','clip3p'], inplace=True)
    df.rename(columns=dict(clip_sta='sta',clip_end='end'),inplace=True)
    df = df[FINAL_COLORDER+['orig']]
    df.reset_index(drop=True, inplace=True)

    #print(df.head(50))
    return df

if __name__ == "__main__":
    args = get_parser().parse_args()
    df = sam2softclippedbed(args.sambam_file)
    df.to_csv(sys.stdout,index=False,header=False,sep='\t')

