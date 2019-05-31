awk -v RS='>' -v FS='\n' '$1~/^C/ {print $0}' raw_data/symbio_phylo/SymPortalDBSequences.fasta > output/symbio_phylo/SymPortalDBSequences_clean.fasta

cat output/symbio_phylo/SymPortalDBSequences_clean.fasta raw_data/symbio_phylo/voolstra_symbio_processed.fasta > output/symbio_phylo/symportal_plus_voolstraDB.fasta

mafft --adjustdirectionaccurately --maxiterate 1000 --ep 0 --genafpair output/symbio_phylo/symportal_plus_voolstraDB.fasta | sed '/>/s/$/</' | tr -d '\n' | sed 's/>/\
>/g' | sed 's/</\
/g' | sed '/^$/d' > output/symbio_phylo/symportal_plus_voolstraDB_aligned.fasta

FastTree -gtr -nt output/symbio_phylo/symportal_plus_voolstraDB_aligned.fasta > output/symbio_phylo/symportal_plus_voolstraDB_aligned_FastTreeGTR.tree
