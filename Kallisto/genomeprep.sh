cat ../fastq/Mus_musculus.GRCm38.cdna.all.fa ../fastq/Mus_musculus_c3hhej.C3H_HeJ_v1.cdna.all.fa ../fastq/Mus_musculus_dba2j.DBA_2J_v1.cdna.all.fa > transcripts2.fa
kallisto index -i transcripts2.idx transcripts2.fa
