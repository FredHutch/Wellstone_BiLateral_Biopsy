#!/bin/bash
#./sort_index.sh
#SBATCH -N1 -n4 -t 4-0 -p campus-new --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

ml purge
ml SAMtools/1.10-GCCcore-8.3.0

bam_file=$1
sample="$(basename -- $bam_file)"
bam_dir="/fh/scratch/delete90/tapscott_s/hg38.BiLat.FSHD.biopsy/sort_index_bam"
sample_name=$(echo "${sample}" | sed 's/trimmed_//' | sed 's/.fastq.gz.subread//')

cd $bam_dir
# sort and index bam
samtools sort -@ 4 $bam_file -o $bam_dir/$sample_name 
samtools index -@ 4 $bam_dir/$sample_name
touch $bam_dir/$sample_name.sort-index-Done.txt
exit 0