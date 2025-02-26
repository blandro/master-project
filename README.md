# master-project
This project aims to improve the annotation of metabolite transporters. Primarily through TCDB, but also with resources as Rhea, ChEBI and UniProt. The goal is to connect genomes to specific properties, like reactions, substrates, rates, locations and more.

## Project Description
The project is split in two methodologies in order to better detect flaws, and analyze strengths and weaknesses of the mappings between databases. These are Approach 1 and Approach 2.

#### Approach 1 
A1 is based solely on what TCDB provides of data. Its data is connected to UniProt, which in turn is mapped to Rhea. The Rhea section acts as a "fact checker" of the mechanisms and substrates that have been extracted from TCDB. Rhea is in fact seen as the ground truth.

#### Approach 2
A2 looks at the same problem, but from the other way around. Starting with all the data on Rhea, this maps to all the assigned accesions in UniProt, before what is available, is mapped to TCDB, and its data, meaning the AA and substrates.

#### Pipelines
Pipeline 1 and Pipeline 2 are just the execution of a BLAST against the databases that was created from A1 and A2. Follow the instructions as provided. After the BLAST, the results (with SubjectID) is merged with a the data the database was built on. I.e. P1 merges the results with data obtained from A1, and vv. for P2. This connects the identified proteins with certain reactions and substrates.

### ChEBI
ChEBI is a dictionary of molecular entities focused on 'small' chemical compounds. All transporters that have any substrate mapped to itself, has a ChEBI ID. However, there exist old instances of ChEBI IDs. Some of these have unfortunately been employed by TCDB. These IDs are called 'secondary IDs'. The ChEBI-folder serves two main purposes. It converts any secondary ID to its primary ID. This is possible as ChEBI has a complete mapping between primary and secondary IDs, due to the arbitrary usage of them, on different sites and softwares. Like TCDB. The second purpose, is the hierarchy creation of ChEBI IDs. If possible for a later part of the project, the amount of desscendants can be used to filter out the non-specific instances of a dataset. A non-specific substrate is 'a molecule'. This has several thousands of descendants. A more specific substrate is 'hydron' (hydrogen proton). This has only three descendants, 'proton', 'deuteron', and 'triton'. However, hydron is specific enough to have data regarding its properties, attributed.

### UniProt
The UniProt-folder contains the relationship between UIDs, RSIDs and RIDs. Obtaining these relations is a cumbersome and large process to do locally. Hence, the queries used to do this, have been provided for the user.

### Misc
This folder is just various quantitative and qualitative analysis of the composition of TCDB, used to determine various aspects of the A1 and A2.