# master-project
This project aims to evaluate how reliably current tools and databases can be used to identify and characterize transporters and their transport reactions. Primarily through TCDB, but also heavily supported and supplemented with Rhea, ChEBI and UniProt. The goal is to connect genomes to specific properties, like reactions, substrates and ChEBI IDs. This is gathered in Benchmarked Data on Curated Transporters (BDCT). Subsequently, BDCT is used to build a BLASTp pipeline, used for a case study of iML1515.

## Project Description
The project is split in two methodologies in order to better detect flaws, and analyze strengths and weaknesses of the mappings between databases. These are Approach 1 and Approach 2, and represent DB1 and DB2 in the thesis. Further on, these were concatinated and polished, giving BDCT (transporters_df.tsv). This is the reference-database used in the pipeline (where it is called transporters_db). Further on, this pipeline was used to compare findings from the genome on which iML1515 is built upon. This was done for multiple combinations of E-values and sequence identies, before the findings were analyzed. Analysis was performed on the results from the case study, but also on BDCT itself.

#### Approach 1 
A1 is based on what TCDB provides of data. Its data is connected to UniProt, which in turn is mapped to Rhea. The Rhea section acts as a "fact checker" of the mechanisms and substrates that have been extracted from TCDB. Rhea is in fact seen as the ground truth. The end product is in essence DB1.

#### Approach 2
A2 looks at the same problem, but from the other way around. Starting with all the data on Rhea, this maps to all the assigned accesions in UniProt, before what is available, is mapped to TCDB, and its data, meaning the AA and substrates. The end product is in essence DB2.

#### 1+2
1+2 is the merging of DB1 and DB2 from A1 and A2, respectively. The database BDCT is the file called transporters_df.tsv, as this is used as the FASTA file to create the reference DB for the BLASTp pipeline. The file is called transporters_df becuase it was made long before the name of BDCT saw the light of day, and changing this, would result in many necessary edits to make scripts and notebooks in different folders compile.

#### Analysis
This folder consists of some analysis of BDCT (transporters_dt.tsv). 

#### Pipeline
This folder contains a set-up file for the CLI tool, as well as the genome which was used for iML1515, to compare transporters and transport reactions.

#### Case Study
The case study of iML1515 is conducted. The pipeline is used, and all results stored in outp. A script was written to automatically analyze all results, to achieve PR-plots on identified genes and transporters, and individual plots of transporter genes, and transport reactions common between the results from iML1515 and the BLASTp. That is found in figs. variants contain some temporary testing and analysis of the optimal BLAST configuration. 

#### Misc
This folder is just various quantitative and qualitative analysis of the composition of TCDB, used to determine various aspects of the A1 and A2. It also contains the manual extractions of reaction mechanisms for 50 families. This file is Misc/All_comp/families_mechanisms_entity_final.tsv.

#### ChEBI
ChEBI contains the necessary info to perform s2p and p2n transformations, as well as an analysis on ChEBI IDs in TCDB. It also provides a hierarchy of all descendants of ChEBI IDs, and a short usage example. A conversion from ChEBI ID to ChEBI name is also given, as it is FAIR to make the reactions readable.

#### UniProt
The UniProt-folder contains the relationship between UIDs, RSIDs and RIDs. Obtaining these relations is a cumbersome and large process to do locally. Hence, the queries used to do this, have been provided for the user.

#### BLAST
This is a draft folder, used in creation of the pipeline, to understand the outputs and formatting of BLAST.

#### test_pipelines
The folder is not of interest or use, it only functioned as an intermediary step. Pipeline 1 and Pipeline 2 are just the execution of a BLAST against the databases that were created from A1 and A2. After the BLAST, the results (with SubjectID) are merged with the data the database was built on. I.e. P1 merges the results with data obtained from A1, and vv. for P2. This connects the identified proteins with certain reactions and substrates.

#### files
This folder contains some files that were though to be used across different folders quite a lot. These are not tracked, as they are large. They can be obtained by following the instructions given in the notebooks where they first are needed.