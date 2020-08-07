The file "abt.rda" contains the dataset in Dethlefsen, Les, and David A. Relman. "Incomplete recovery and individualized responses of the human distal gut microbiota to repeated antibiotic perturbation." Proceedings of the National Academy of Sciences 108.Supplement 1 (2011): 4554-4561.

The data is in the form of a phyloseq object (see the "phyloseq" R package).

The file "preprocess.R" contains the code to get the datasets we used in Section 3.2 from the original dataset "abt.rda". The preprocessed dataset contains 3 OTU tables for the 3 subjects in the study and a trimmed phylogenetic tree. The OTUs are aggregated to level 5 in the taxonomic table.
