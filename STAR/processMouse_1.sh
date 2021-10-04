mkdir M1
cd M1
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../../fastq/1-15_S15_R1_001_val_1-cA10-p.fastq ../../fastq/1-15_S15_R2_001_val_2-cT10-p.fastq --quantMode GeneCounts
cd ..

mkdir M2
cd M2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../../fastq/1-16_S16_R1_001_val_1-cA10-p.fastq ../../fastq/1-16_S16_R2_001_val_2-cT10-p.fastq --quantMode GeneCounts
cd ..

mkdir M3
cd M3
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../../fastq/1-17_S17_R1_001_val_1-cA10-p.fastq ../../fastq/1-17_S17_R2_001_val_2-cT10-p.fastq --quantMode GeneCounts
cd ..

mkdir M4
cd M4
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../../fastq/1-18_S18_R1_001_val_1-cA10-p.fastq ../../fastq/1-18_S18_R2_001_val_2-cT10-p.fastq --quantMode GeneCounts
cd ..

mkdir M5
cd M5
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../../fastq/1-7_S13_R1_001_val_1-cA10-p.fastq ../../fastq/1-7_S13_R2_001_val_2-cT10-p.fastq --quantMode GeneCounts
cd ..

mkdir M6
cd M6
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../../fastq/1-8_S14_R1_001_val_1-cA10-p.fastq ../../fastq/1-8_S14_R2_001_val_2-cT10-p.fastq --quantMode GeneCounts
cd ..

mkdir O1
cd O1
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../M1/Unmapped.out.mate1 ../M1/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir O2
cd O2
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../M2/Unmapped.out.mate1 ../M2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir O3
cd O3
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../M3/Unmapped.out.mate1 ../M3/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir O4
cd O4
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../M4/Unmapped.out.mate1 ../M4/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir O5
cd O5
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../M5/Unmapped.out.mate1 ../M5/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir O6
cd O6
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../M6/Unmapped.out.mate1 ../M6/Unmapped.out.mate2 --quantMode GeneCounts
cd ..
