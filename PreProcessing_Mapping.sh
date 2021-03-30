#!/PreProcessing_Mapping.sh

#The present Bash script run the trimming step and check the quality of the reads before and after the trimming phase. 
#Then, it runs the STAR tool to align the reads on the species' genome. 


###1) Pre-processing analysis

#Raw FastQ Data Quality Control
$	fastqc 0_Raw_Data/*fastq -o 0_Raw_Data/QC
$	multiqc 0_Raw_Data/QC -o 0_Raw_Data/QC

#Raw FastQ Data Trimming
$	for f in 0_Raw_Data/*fastq; do trimmomatic SE $f 1_Trimmed_Data/${f%%fastq}"trim.fastq" ILLUMINACLIP:Tool/Adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50; done

#Trimmed FastQ Data Quality Control
$	fastqc 1_Trimmed_Data/*fastq -o 1_Trimmed_Data/QC
$	multiqc 1_Trimmed_Data/QC -o 1_Trimmed_Data/QC


###6) Mapping & Quantification analysis

#Homo sapiens
$	for f in 1_Trimmed_Data/*fastq; do star --genomeDir Reference/Hsapiens/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --outFileNamePrefix 2_Hsapiens_Mapped/${f%%fastq}; done
$	mv 2_Hsapiens_Mapped/*UnMapped* 2_Hsapiens_UnMapped/
$	mv 2_Hsapiens_Mapped/*ReadsPerGene* 3_Quantification/Hsapiens

#Streptococcus pyogenes AP1
$	for f in 2_Hsapiens_UnMapped/*; do star --genomeDir Reference/Spyogenes_AP1/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix 2_Spyogenes_AP1_Mapped/${f%%_Unmapped.out.mate1}; done
$	mv 2_Spyogenes_AP1_Mapped/*ReadsPerGene* 3_Quantification/Spyogenes_AP1

#Streptococcus pyogenes NZ131
$	for f in 2_Hsapiens_UnMapped/*mate1; do star --genomeDir Reference/Spyogenes_NZ131/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix 2_Spyogenes_NZ131_Mapped/${f%%_Unmapped.out.mate1}; done
$	mv 2_Spyogenes_NZ131_Mapped/*ReadsPerGene* 3_Quantification/Spyogenes_NZ131

#Influenza A virus
$	for f in 2_Hsapiens_UnMapped/*mate1; do star --genomeDir Reference/IAV/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix 2_IAV_Mapped/${f%%_Unmapped.out.mate1}; done
$	mv 2_IAV_Mapped/*ReadsPerGene* 3_Quantification/IAV
