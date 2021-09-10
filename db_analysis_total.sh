#!/bin/bash

set -eu -o pipefail

################################################################################
### Arg parse and check                                                      ###
################################################################################
unset db_name
unset local_db_path
unset indir
unset outdir
unset n_threads


#Assign arg variables
while getopts n:l:i:o:t: opt;
do
    case "${opt}" in
        n) db_name=${OPTARG};;
        l) local_db_path=${OPTARG};;
        i) indir=${OPTARG};;
        o) outdir=${OPTARG};;
        t) n_threads=${OPTARG:-1};;
        *)
                        echo 'Unrecognized arguments.' >&2
                        exit 1
    esac
done

#Check that mandatory args are assigned
shift "$(( OPTIND - 1))"
: "${db_name:?Required set db name with -n}"
: "${local_db_path:?Required set local db path with -l}"
: "${indir:?Required set input directory with -i}"
: "${outdir:?Required set output directory with -o}"

################################################################################


### prepare directories
mkdir -p "$outdir"
cat /dev/null > ./commands.txt

#creates for each assembly a folder with the same name as the assembly file
for genome in "$indir"/*.fasta;
do
   mkdir ${genome%.*}; mv "$genome" ${genome%.*};
done;

#for each assembly in the folder finds the ORFs using Prodigal
# (-i, --input file, -d, output --mRNA file) and blastx them against the VFDB
# proteins database using diamond (the database from VFDB is protein-based
# therefore blastx is used, --evalue specifies the threshold for discarding
#  sequences, --outfmt specifies the format, in this case as a table, choice 6,
# including all the options following the command, for every sequence a maximum
# of 25 hits is set as specified by --max-target-seqs
for genome in "$indir"/*/*.fasta;
do
  echo -e "prodigal -i "$genome" -d "${genome%.*}"_orfs.fna; diamond blastx --query "${genome%.*}"_orfs.fna --db $local_db_path --evalue 1e-06 --outfmt 6 qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore full_qseq --max-target-seqs 25 --out "${genome%.*}"_"$db_name"_hits.out6" >> commands.txt;
done;

#runs the process in parallel
parallel -j "$n_threads" < ./commands.txt

# Delete DS_Store files from mac and windows file systems
find "$indir" -name '.DS_Store' -type f -delete

#Python script to sort the hits by percentage of identity, filter out those
#below 80% similarity, and compile a single master table with all the hits for
#all the assembly
python filter_sort_concat_vfdb.py -i "$indir" -db "$db_name" -o "$outdir"

awk -F , 'NR>1{print ">"$17"\n"$16}' "$outdir"/total_"$db_name"_table.csv  > "$outdir"/total_"$db_name"_table.fasta

mkdir -p "$outdir"/"$db_name"_data

vsearch --cluster_size "$outdir"/total_"$db_name"_table.fasta --notrunclabels --clusters "$outdir"/"$db_name"_data/"$db_name"Cluster --id 0.95

for i in "$outdir"/"$db_name"_data/*;
do
  mv $i $(echo $i".fasta");
done ;

python gloome_to_csv.py -g "$outdir"/"$db_name"_data/ -o "$outdir"/output_presence_absence

#To run gainLoss, uncomment
#gainLoss paramfile.txt
