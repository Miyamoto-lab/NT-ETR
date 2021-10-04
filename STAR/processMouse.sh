
#First run mouse analysis
cd basic_aAm1_S12_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm1_R1_trimmed.fastq ../aAm1_R2_trimmed.fastq --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../mousey.gtf > Res.counts
cd ..

cd basic_aAm2_S13_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm2_R1_trimmed.fastq ../aAm2_R2_trimmed.fastq --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../mousey.gtf > Res.counts
cd ..

cd basic_aAm3_S14_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm3_R1_trimmed.fastq ../aAm3_R2_trimmed.fastq --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../mousey.gtf > Res.counts
cd ..

cd basic_1_S9_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT1_R1_trimmed.fastq ../NT1_R2_trimmed.fastq --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../mousey.gtf > Res.counts
cd ..

cd basic_2_S10_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT2_R1_trimmed.fastq ../NT2_R2_trimmed.fastq --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../mousey.gtf > Res.counts
cd ..

cd basic_3_S11_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT3_R1_trimmed.fastq ../NT3_R2_trimmed.fastq --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../mousey.gtf > Res.counts
\cd ..

#Now Oryx
cd basic_1_S9_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_1_S9_2/Unmapped.out.mate1 ../basic_1_S9_2/Unmapped.out.mate2 --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../oryx.gtf > Res.counts
cd ..

cd basic_2_S10_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_2_S10_2/Unmapped.out.mate1 ../basic_2_S10_2/Unmapped.out.mate2 --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../oryx.gtf > Res.counts
cd ..

cd basic_3_S11_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_3_S11_2/Unmapped.out.mate1 ../basic_3_S11_2/Unmapped.out.mate2 --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../oryx.gtf > Res.counts
cd ..

cd basic_aAm1_S12_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm1_S12_2/Unmapped.out.mate1 ../basic_aAm1_S12_2/Unmapped.out.mate2 --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../oryx.gtf > Res.counts
cd ..

cd basic_aAm2_S13_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm2_S13_2/Unmapped.out.mate1 ../basic_aAm2_S13_2/Unmapped.out.mate2 --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../oryx.gtf > Res.counts
cd ..

cd basic_aAm3_S14_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm3_S14_2/Unmapped.out.mate1 ../basic_aAm3_S14_2/Unmapped.out.mate2 --quantMode GeneCounts
htseq-count -m intersection-nonempty -i gene_id Aligned.out.sam ../oryx.gtf > Res.counts
cd ..
