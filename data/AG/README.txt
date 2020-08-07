This folder contains the American Gut data for the case studies. The two datasets "AG_ibd.RData" and "AG_diabetes.RData" are preprocessed based on the July 29, 2016 version of the fecal data from the American Gut Project. Each dataset contains an OTU table and a phylogenetic tree. The OTU tables contain counts of the top 75 OTUs selected by the total counts of the OTUs in all samples. 

AG_ibd.RData is the IBD dataset used in Section 4.1. The dataset has 189 samples from IBD patients in the study.

AG_diabetes.RData is the Diabetes dataset used in Section 4.2. The dataset has 106 samples from diabetes patients in the study.

These two datasets are obtained by directly running the preprocessing code in "Data_process_ibd.R" and "Data_process_diabetes.R", respectively. The code preprocess the original AG data by selecting subset of samples and OTUs. Details of the preprocessing procedure is given in the paper.

The original July 29, 2016 version of the fecal dataset in the American Gut Project contains 4 parts:
1. "97_otu_taxonomy.txt" contains a taxa table of all OTUs in the study;
2. "97_otus.tree" contains a rooted full binary phylogenetic tree of these OTUs;
3. "ag_fecal_from_biom.txt" contains the entire OTU table of all samples. This table is directly obtained based on the sequencing reads.
4. "ag_fecal.txt" is the metadata of the study.



