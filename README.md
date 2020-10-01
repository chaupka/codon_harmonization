# codon_harmonization
Create/Convert and visualize codon usage tables based on genbank files ('gb') or existent tables from 'kazusa.or.jp'. For generation of codon usage tables via genbank, consider only CDS's with at least 300 bp. Optimize the codon usage of your gene not by maximizing the codon frequency, but by adapting a combination of similar codon frequency and relative adaptiveness values present in the original host organism. Change startcodons to 'ATG'.

# Notes:
For conversion of codon usage tables from kazusa.or.jp you need to check the 'A style like CodonFrequency output in GCG Wisconsin PackageTM' option, then copy and save the output as a 'csv' file.
Requires Python3 and modules:
- Biopython
- matplotlib