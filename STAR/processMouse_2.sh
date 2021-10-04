mkdir basic_aAm1_S12_2B
cd basic_aAm1_S12_2B
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm1_R1_trimmed.fastq ../aAm1_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_aAm2_S13_2B
cd basic_aAm2_S13_2B
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm2_R1_trimmed.fastq ../aAm2_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_aAm3_S14_2B
cd basic_aAm3_S14_2B
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm3_R1_trimmed.fastq ../aAm3_R2_trimmed.fastq --quantMode GeneCounts
cd ..


mkdir basic_1_S9_2B
cd basic_1_S9_2B
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT1_R1_trimmed.fastq ../NT1_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_2_S10_2B
cd basic_2_S10_2B
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT2_R1_trimmed.fastq ../NT2_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_3_S11_2B
cd basic_3_S11_2B 
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT3_R1_trimmed.fastq ../NT3_R2_trimmed.fastq --quantMode GeneCounts
cd ..


mkdir basic_aAm1_S12_2
cd basic_aAm1_S12_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm1_R1_trimmed.fastq ../aAm1_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_aAm2_S13_2
cd basic_aAm2_S13_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm2_R1_trimmed.fastq ../aAm2_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_aAm3_S14_2
cd basic_aAm3_S14_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm3_R1_trimmed.fastq ../aAm3_R2_trimmed.fastq --quantMode GeneCounts
cd ..


mkdir basic_1_S9_2
cd basic_1_S9_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT1_R1_trimmed.fastq ../NT1_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_2_S10_2
cd basic_2_S10_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT2_R1_trimmed.fastq ../NT2_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_3_S11_2
cd basic_3_S11_2
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT3_R1_trimmed.fastq ../NT3_R2_trimmed.fastq --quantMode GeneCounts
cd ..

#Oryx ... nnt done
mkdir basic_1_S9_2C
cd basic_1_S9_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_1_S9_2/Unmapped.out.mate1 ../basic_1_S9_2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_2_S10_2C
cd basic_2_S10_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_2_S10_2/Unmapped.out.mate1 ../basic_2_S10_2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_3_S11_2C
cd basic_3_S11_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_3_S11_2/Unmapped.out.mate1 ../basic_3_S11_2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm1_S12_2C
cd basic_aAm1_S12_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm1_S12_2/Unmapped.out.mate1 ../basic_aAm1_S12_2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm2_S13_2C
cd basic_aAm2_S13_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm2_S13_2/Unmapped.out.mate1 ../basic_aAm2_S13_2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm3_S14_2C
cd basic_aAm3_S14_2C
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm3_S14_2/Unmapped.out.mate1 ../basic_aAm3_S14_2/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

#Oryx
mkdir basic_1_S9_2D
cd basic_1_S9_2D
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_1_S9_2B/Unmapped.out.mate1 ../basic_1_S9_2B/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_2_S10_2D
cd basic_2_S10_2D
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_2_S10_2B/Unmapped.out.mate1 ../basic_2_S10_2B/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_3_S11_2D
cd basic_3_S11_2D
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_3_S11_2B/Unmapped.out.mate1 ../basic_3_S11_2B/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm1_S12_2D
cd basic_aAm1_S12_2D
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm1_S12_2B/Unmapped.out.mate1 ../basic_aAm1_S12_2B/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm2_S13_2D
cd basic_aAm2_S13_2D
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm2_S13_2B/Unmapped.out.mate1 ../basic_aAm2_S13_2B/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm3_S14_2D
cd basic_aAm3_S14_2D
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm3_S14_2B/Unmapped.out.mate1 ../basic_aAm3_S14_2B/Unmapped.out.mate2 --quantMode GeneCounts
cd ..



#Base2

mkdir basic_aAm1_S12_2E
cd basic_aAm1_S12_2E
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm1_R1_trimmed.fastq ../aAm1_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_aAm2_S13_2E
cd basic_aAm2_S13_2E
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm2_R1_trimmed.fastq ../aAm2_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_aAm3_S14_2E
cd basic_aAm3_S14_2E
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../aAm3_R1_trimmed.fastq ../aAm3_R2_trimmed.fastq --quantMode GeneCounts
cd ..


mkdir basic_1_S9_2E
cd basic_1_S9_2E
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT1_R1_trimmed.fastq ../NT1_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_2_S10_2E
cd basic_2_S10_2E
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT2_R1_trimmed.fastq ../NT2_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_3_S11_2E
cd basic_3_S11_2E
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40  --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../NT3_R1_trimmed.fastq ../NT3_R2_trimmed.fastq --quantMode GeneCounts
cd ..

mkdir basic_1_S9_2F
cd basic_1_S9_2F
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_1_S9_2E/Unmapped.out.mate1 ../basic_1_S9_2E/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_2_S10_2F
cd basic_2_S10_2F
STAR --genomeDir  ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_2_S10_2E/Unmapped.out.mate1 ../basic_2_S10_2E/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_3_S11_2F
cd basic_3_S11_2F
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_3_S11_2E/Unmapped.out.mate1 ../basic_3_S11_2E/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm1_S12_2F
cd basic_aAm1_S12_2F
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm1_S12_2E/Unmapped.out.mate1 ../basic_aAm1_S12_2E/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm2_S13_2F
cd basic_aAm2_S13_2F
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm2_S13_2E/Unmapped.out.mate1 ../basic_aAm2_S13_2E/Unmapped.out.mate2 --quantMode GeneCounts
cd ..

mkdir basic_aAm3_S14_2F
cd basic_aAm3_S14_2F
STAR --genomeDir ../../Oryx/MGenome --sjdbGTFfile ../mousey.gtf --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm3_S14_2E/Unmapped.out.mate1 ../basic_aAm3_S14_2E/Unmapped.out.mate2 --quantMode GeneCounts
cd ..


#Oryx ... nnt done
mkdir basic_1_S9_2G
cd basic_1_S9_2G
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_1_S9_2/Unmapped.out.mate1 ../basic_1_S9_2/Unmapped.out.mate2 --qu$
cd ..

mkdir basic_2_S10_2G
cd basic_2_S10_2G
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_2_S10_2/Unmapped.out.mate1 ../basic_2_S10_2/Unmapped.out.mate2 --$
cd ..

mkdir basic_3_S11_2G
cd basic_3_S11_2G
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_3_S11_2/Unmapped.out.mate1 ../basic_3_S11_2/Unmapped.out.mate2 --$
cd ..

mkdir basic_aAm1_S12_2G
cd basic_aAm1_S12_2G
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm1_S12_2/Unmapped.out.mate1 ../basic_aAm1_S12_2/Unmapped.out.ma$
cd ..

mkdir basic_aAm2_S13_2G
cd basic_aAm2_S13_2G
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm2_S13_2/Unmapped.out.mate1 ../basic_aAm2_S13_2/Unmapped.out.ma$
cd ..

mkdir basic_aAm3_S14_2G
cd basic_aAm3_S14_2G
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm3_S14_2/Unmapped.out.mate1 ../basic_aAm3_S14_2/Unmapped.out.ma$
cd ..

mkdir basic_1_S9_2H
cd basic_1_S9_2H
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_1_S9_2B/Unmapped.out.mate1 ../basic_1_S9_2B/Unmapped.out.mate2 --$
cd ..

mkdir basic_2_S10_2H
cd basic_2_S10_2H
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_2_S10_2B/Unmapped.out.mate1 ../basic_2_S10_2B/Unmapped.out.mate2 $
cd ..

mkdir basic_3_S11_2H
cd basic_3_S11_2H
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_3_S11_2B/Unmapped.out.mate1 ../basic_3_S11_2B/Unmapped.out.mate2 $
cd ..

mkdir basic_aAm1_S12_2H
cd basic_aAm1_S12_2H
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm1_S12_2B/Unmapped.out.mate1 ../basic_aAm1_S12_2B/Unmapped.out.$
cd ..

mkdir basic_aAm2_S13_2H
cd basic_aAm2_S13_2H
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm2_S13_2B/Unmapped.out.mate1 ../basic_aAm2_S13_2B/Unmapped.out.$
cd ..

mkdir basic_aAm3_S14_2H
cd basic_aAm3_S14_2H
STAR --genomeDir ../../Oryx/OGeonome2 --sjdbGTFfile ../oryx.gtf --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --outReadsUnmapped Fastx --sjdbOverhang 149 --readFilesIn ../basic_aAm3_S14_2B/Unmapped.out.mate1 ../basic_aAm3_S14_2B/Unmapped.out.$
cd ..
