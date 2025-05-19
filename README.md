# master-project
This project aims to improve the annotation of metabolite transporters. Primarily through TCDB, but also with resources as Rhea, ChEBI and UniProt. The goal is to connect genomes to specific properties, like reactions, substrates and ChEBI IDs. This is gathered in Transporter Database (TDB). Subsequently, TDB is used to build a BLASTp pipeline, used for a case study of iML1515.

## Project Description
The project is split in two methodologies in order to better detect flaws, and analyze strengths and weaknesses of the mappings between databases. These are Approach 1 and Approach 2.

#### Approach 1 
A1 is based solely on what TCDB provides of data. Its data is connected to UniProt, which in turn is mapped to Rhea. The Rhea section acts as a "fact checker" of the mechanisms and substrates that have been extracted from TCDB. Rhea is in fact seen as the ground truth. The end product is in essence DB1.

#### Approach 2
A2 looks at the same problem, but from the other way around. Starting with all the data on Rhea, this maps to all the assigned accesions in UniProt, before what is available, is mapped to TCDB, and its data, meaning the AA and substrates. The end product is in essence DB2.

#### 1+2
1+2 is the merging of DB1 and DB2 from A1 and A2, respectively. The database is the file called transporters_df.tsv, as this is used as the FASTA file to create the reference DB for the BLASTp pipeline.

#### Pipeline
This folder contains a set-up file for the CLI tool, as well as the proteome which was used for iML1515, to compare transporters and transport reactions.

#### test_pipelines
Pipeline 1 and Pipeline 2 are just the execution of a BLAST against the databases that was created from A1 and A2. Follow the instructions as provided. After the BLAST, the results (with SubjectID) is merged with a the data the database was built on. I.e. P1 merges the results with data obtained from A1, and vv. for P2. This connects the identified proteins with certain reactions and substrates.

#### Case Study
The case study of iML1515 is conducted. The pipeline is used, and all results stored in outp. A script was written to automatically analyze all results, to achieve both a PR-plot, and individual plots of transporter genes, and transport reactions common between the results from iML1515 and the BLASTp. That is found in figs. variants contain some temporary testing. 

#### BLAST
This is just a draft folder, used in creation of the pipeline.

### ChEBI
ChEBI is a dictionary of molecular entities focused on 'small' chemical compounds. All transporters that have any substrate mapped to itself, has a ChEBI ID. However, there exist old instances of ChEBI IDs. Some of these have unfortunately been employed by TCDB. These IDs are called 'secondary IDs'. The ChEBI-folder serves two main purposes. It converts any secondary ID to its primary ID. This is possible as ChEBI has a complete mapping between primary and secondary IDs, due to the arbitrary usage of them, on different sites and softwares. Like TCDB. The second purpose, is the hierarchy creation of ChEBI IDs. If possible for a later part of the project, the amount of desscendants can be used to filter out the non-specific instances of a dataset. A non-specific substrate is 'a molecule'. This has several thousands of descendants. A more specific substrate is 'hydron' (hydrogen proton). This has only three descendants, 'proton', 'deuteron', and 'triton'. However, hydron is specific enough to have data regarding its properties, attributed.

### UniProt
The UniProt-folder contains the relationship between UIDs, RSIDs and RIDs. Obtaining these relations is a cumbersome and large process to do locally. Hence, the queries used to do this, have been provided for the user.

### Misc
This folder is just various quantitative and qualitative analysis of the composition of TCDB, used to determine various aspects of the A1 and A2. It also contains the manual extractions of reaction mechanisms for 50 families.