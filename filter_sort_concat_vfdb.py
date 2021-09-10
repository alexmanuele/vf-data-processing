import pandas as pd
import glob
from pathlib import Path
import argparse
import sys

#Suppose the database is defined
# TODO make this an arg
#TODO proper argparse and usage check
argparser = argparse.ArgumentParser(description='Process results from diamond.')
argparser.add_argument('-i', help='Input directory containing relevant files/directories')
argparser.add_argument('-db', help='Name of database, e.g. vfdb or bacmet.')
argparser.add_argument('-o', help='Ouput directory.')

if __name__ == "__main__":

    args = argparser.parse_args()

    db = args.db
    inpath = args.i
    outpath = args.o

    all_VFDB_results = []

    for file in Path(inpath).glob("*"):
        stem = str(file).replace(inpath, '').replace('_assembly', '')
        hit_dir = str(file) + '/' + stem + '_{0}_hits.out6'.format(db)

        df = pd.read_csv(hit_dir, sep='\t', names=["qseqid", "sseqid", "pident", "slen", "qlen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sequence"]) #reads each output file in every folder and creates a table with those titles
        df['Genome'] = stem
        df['qcover'] = df['length']/df['slen']
        unique_hits = df.sort_values('pident', ascending=False).drop_duplicates(['qseqid', 'sequence'])
        filtered_unique_hits = unique_hits[(unique_hits['pident'] >= 60) & (unique_hits['qcover'] >= 0.6)]
        hit_dir2 = "{0}/{1}".format(file, 'filtered_unqiue_hits2.csv')
        filtered_unique_hits.to_csv(hit_dir2)
        all_VFDB_results.append(filtered_unique_hits)

    all_VFDB_results = pd.concat(all_VFDB_results)

    all_VFDB_results.to_csv(outpath + '/total_VFDB_table.csv')
