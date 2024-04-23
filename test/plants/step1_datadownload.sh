# Download the genome sequence files, genome length files, TE annotation files, and genome annotation files required to run teena.
# Download: from Phytozome.
# This is an example for data download of rice(oryza_sativa:IRGSP1).
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.58.gff3.gz

gunzip *gz
