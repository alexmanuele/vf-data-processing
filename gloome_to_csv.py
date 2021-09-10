from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict


def check_folder(path):
    """
    check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_dir():
        return path
    else:
        raise argparse.argumenttypeerror("{path} isn't a folder")


def build_presence_absence_string(gene_folder):
    """
    parses gene fasta files and creates a presence absence .csv files
    where each genome is represented by the accession
    and each gene by either a 1 or 0 depending on presence or absence in
    that genome
    """
    # check if gene fasta files can be found
    gene_fasta_files = list(gene_folder.glob("*.fasta"))
    if len(gene_fasta_files) == 0:
        raise valueerror("no .fasta files in {gene_folder}")

    # * * * 
    # instead of printing the gene order, we will keep a list. This will become
    # the columns of our csv.
    gene_order = [gene.stem for gene in gene_fasta_files]

    # parse gene fasta and create a genome -> gene dictionary
    genome_to_gene_dict = defaultdict(list)
    for gene_fasta in gene_fasta_files:
        gene_name = gene_fasta.stem
        for record in SeqIO.parse(str(gene_fasta), 'fasta'):
            genome_to_gene_dict[record.id].append(gene_name)

    # * * *
    # Instead of a fasta-style alignment file, we will build .csvs
    # We'll make a list of records, and build the records into a csv.
    wide_records = []
    long_records = []

    for genome, genes in genome_to_gene_dict.items():
        presence_absence_string = ""
        # * * * Keep track of all the genes per genome,.
        wide_record = []
        for gene in gene_fasta_files:
            # * * *
            # Instead of string, use int and build the csv.
            present = 0
            if gene.stem in genes:
                present = 1

            long_records.append({'genome': genome, 'gene': gene.stem, 'presence':present})
            wide_record.append((gene.stem, present))
        # * * *
        # Build the wide record
        wide = {'genome':genome}
        for tup in wide_record:
            wide[tup[0]] = tup[1]
        wide_records.append(wide)


    wide_df = pd.DataFrame.from_records(wide_records)
    long_df = pd.DataFrame.from_records(long_records)
    return wide_df, long_df

if __name__ == '__main__':

    parser = argparse.ArgumentParser("generate gain/loss alignment from a "
                                     "folder of per-gene fasta files")
    parser.add_argument("-g", "--gene_folder", required=True,
                        type=check_folder,
                        help="folder containing *.fasta for each gene with "
                              "genome indicated by accession")

    # * * *
    # Remove the file extension from the output argument.
    parser.add_argument("-o", "--output",
                        default="output_presence_absence",
                        help="output path for presence absence alignment")

    args = parser.parse_args()

    wide_df, long_df = build_presence_absence_string(args.gene_folder)

    print("writing alignment to {args.output}")
    # * * * 
    # Use a format string so we can write both csv files.
    # We use pandas to handle the file formatting.
    wide_df.to_csv("{}_wide.csv".format(args.output), sep=',', index=False)
    long_df.to_csv("{}_long.csv".format(args.output), sep=',', index=False)
