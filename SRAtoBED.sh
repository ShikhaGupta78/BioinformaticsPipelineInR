# Creation of FastQ file from SRA file, first parameter is the path of the fastq-dump and the second parameter is the location of the SRA file
/home/tools/sratoolkit.2.5.1-centos_linux64/bin/fastq-dump.2.5.1 /PathToSRA/FolderName/SRA/$1

# BWA command to map the generated FastQ file to the reference genome, and index SAI file is generated as the output
bwa aln -t 8 /media/lpmb3/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa /media/lpmb3/Shikha/SRA/$1.fastq > $1.sai

# BWA command to finally generate the SAM file from the index SAI and FastQ files
bwa samse /media/lpmb3/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa $1.sai $1.fastq > $1.aln.sam

# SAMtools command to convert SAM file format to binary version BAM file format
samtools view -Shu $1.aln.sam > $1.bam

# SAMtools command to sort the generated BAM file
samtools sort $1.bam $1_sorted.bam

# Generate index file out of the sorted BAM file
samtools index $1_sorted.bam $1_sorted.bam.bai

# Command to generate the Q peak files from the sorted BAM file
/home/tools/Q/Q -l 18 -x 9 -t $1_sorted.bam -o $1.peaks

# Command to concatenate all bed files generated from previous Q command into 1 file with output file name ALL_total_peaks.bed
cat *.bed > ALL_total_peaks.bed

# Command to sort the ALL_total_peaks.bed file and output file generated is named ALL_total_peaks_sorted.bed
bedtools sort -i ALL_total_peaks.bed>ALL_total_peaks_sorted.bed

# Command to merge the previously sorted bed file ALL_total_peaks_sorted.bed
bedtools merge -i ALL_total_peaks_sorted.bed > ALL_total_peaks_sorted_merged.bed 

# bedtools multicov command
bedtools multicov -q 20 -bams SRR5445252_sorted.bam SRR5445251_sorted.bam SRR5445213_sorted.bam SRR5445212_sorted.bam -bed /media/lpmb3/Shikha/Q/testShikhaFinalest.bed >
/media/lpmb3/Shikha/FolderName/final_BED.txt & disown

