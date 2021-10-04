samtools faidx oryx.1_HiC.fasta
cut -f1,2 oryx.1_HiC.fasta.fai > oryx.chrom.sizes
grep -f uplist.txt  TSS.bed > upTSS.bed
grep -f downlist.txt TSS.bed > downTSS.bed

awk '{ print $2 "\t" $3 "\t" $4 "\t" $6}' upTSS.bed > upTSSa.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $6}' downTSS.bed > downTSSa.bed


bedtools slop -i upTSSa.bed -g ../Oryx/oryx.chrom.sizes  -b 1000 > upTSSs.bed
bedtools slop -i upTSSa.bed -g ../Oryx/oryx.chrom.sizes  -b 1000 > downTSSs.bed

fastaFromBed -fi ../Oryx/oryx.1_HiC.fasta -bed upTSSs.bed -fo upTSSs.bed.fasta
findMotifs.pl upTSSs.bed.fasta fasta UP -mset vertebrates

fastaFromBed -fi ../Oryx/oryx.1_HiC.fasta -bed downTSSs.bed -fo downTSSs.bed.fasta
findMotifs.pl downTSSs.bed.fasta fasta DOWN -mset vertebrates
