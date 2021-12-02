# Beat AML Proteomics analysis
This repository houses the exploratory work that we are doing to evaluate the role of proteomics and phosphoproteomics in measuring patient diversity in Acute Myeloid Leukemia, both in the patient outcome and the response of patient samples to ~ex vivo~ drugs.


## Repository Overview

This repository contains scripts that pull data from a [Synapse repository](http://synapse.org/ptrc) to carry out the various analysis steps needed. You will need to acquire a [synapse username](http://synapse.org/register) to access synapse, and become a [certified user](https://docs.synapse.org/articles/accounts_certified_users_and_profile_validation.html) to add data, but after that you will be set with future projects. You will then need to navigate to the [PNNL/OHSU Synapse page](http://synapse.org/ptrc) to request access.


### Before you begin...

This repository is only the basic structure of the tools needed, not the end-to-end analysis. Here are the steps you'll need to use this:

1- Read up on the tools
  - GitHub requires some basic protocols such as pull requests and commits, you should try to get a basic understanding. I found this [tutorial](https://medium.com/@jonathanmines/the-ultimate-github-collaboration-guide-df816e98fb67) that can serve as a starting point.
  - Synapse also has a bit of learning curve. To understand what Synapse is1 and isn't, check out [this document](https://docs.synapse.org/articles/getting_started.html).
2- Get [RStudio](http://rstudio.org). Basic R is essential, but RStudio will make your life a lot easier, I promise!
3- Install the [synapse python client](https://python-docs.synapse.org/build/html/index.html), create a [`.synapseConfig` file](https://python-docs.synapse.org/build/html/Credentials.html) in your home directory.
4- Clone this repository - it has all that you will need to contribute and run this analysis.

## Beat AML Processing
Here we descirbe the processing of the BeatAML Data

### Proteomics and phosphoproteomics processing

This repository contains the code for the normalization and processing.

This is derived from the P3 proteomics workflow and can be found in the [proteomics](./proteomics) folder. It contains the scripts required to process and normalize the data. It also contains the [study design](./proteomics/study_design) files that are required to do the processing.

Once the data is processed from DMS it is uploaded to Synapse in the [Proteomics and Quality Control](https://www.synapse.org/#!Synapse:syn24171150) folder.

The files are all stored in Synapse so that they can be downloaded and shared.
| Description | Link |
| --- | --- |
| Global proteomics data files | [syn25714186](https://www.synapse.org/#!Synapse:syn25714186) |
| Phosphoproteomics data files | [syn25714185](https://www.synapse.org/#!Synapse:syn25714185) |
| Metadata file | [syn25807733](https://www.synapse.org/#!Synapse:syn25807733) |

The data was pushed from raw files to long-form tables for facile querying and viewing:

| Description | Normalization/filtering| Link |
| --- | --- | --- |
| Global Proteomics | | |
| Global proteomics |  | |
| Global phosphoproteomics | Camilo add here |[syn26477193](https://www.synapse.org/#!Synapse:syn26477193) |
| Global phosphoproteomics |||
| Global phosphoprotoemics |||
| Global phosphoproteomics |||

### Gene mutation, RNASeq, Clinical data

These data are also uploaded to Synapse, and then parsed.

| Description | Link |
| --- | -- |
| Waves 1to4 WES Data | [syn2648827](https://www.synapse.org/#!Synapse:syn26428827/tables/) |
| Waves 1to4 RNASeq Data | [syn26428813](https://www.synapse.org/#!Synapse:syn26428813) |


## Data Analysis

This section aspirationally aims to serve as the outline for the manuscript we are building.

### Cohort exploration and summarization
The first figure of the manuscript will require visualizing the Beat AML cohort and the data we have.


### Mutational profiling

### Immune deconvolution

### Functional analysis
