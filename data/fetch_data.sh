curl -O https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz && \
    gunzip -c Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    > TAIR10.fasta

curl -O https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gff3.gz && \
    gunzip -c Arabidopsis_thaliana.TAIR10.56.gff3.gz \
    > TAIR10.gff3



rm Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
rm Arabidopsis_thaliana.TAIR10.56.gff3.gz
