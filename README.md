# NEO2IS_ITHS
# A novel integrated approach to predicting cancer immunotherapy efficacy
# Codes
### Somatic mutations and structural variants were detected using WES and RNA-seq data. Translated proteins were chopped into 8-11 kmers peptides encompassing the mutated residues (for SNVs and indels) and fusion breakpoint (for fusions) until a stop codon. Mutated peptides with binding affinity <0.5 % Rank determined by NetMHCpan were considered as candidate neoantigens. The neoantigen load score was calculated based on outputs of deepHLApan model
