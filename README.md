# Beat AML Proteomics analysis
This repository houses the exploratory work that we are doing to evaluate the role of proteomics and phosphoproteomics in measuring patient diversity in Acute Myeloid Leukemia, both in the patient outcome and the response of patient samples to ~ex vivo~ drugs.


## Repository Overview

This repository contains scripts that pull data from a [Synapse repository](http://synapse.org/ptrc) to carry out the various analysis steps needed. You will need to acquire a [synapse username](http://synapse.org/register) to access synapse, and become a [certified user](https://docs.synapse.org/articles/accounts_certified_users_and_profile_validation.html) to add data, but after that you will be set with future projects. You will then need to navigate to the [PNNL/OHSU Synapse page](http://synapse.org/ptrc) to request access.


### Before you begin...

This repository is only the basic structure of the tools needed, not the end-to-end analysis. Here are the steps you'll need to use this:

1. Read up on the tools.  GitHub requires some basic protocols such as pull requests and commits, you should try to get a basic understanding. I found this [tutorial](https://medium.com/@jonathanmines/the-ultimate-github-collaboration-guide-df816e98fb67) that can serve as a starting point.  Synapse also has a bit of learning curve. To understand what Synapse is1 and isn't, check out [this document](https://docs.synapse.org/articles/getting_started.html).
2. Get [RStudio](http://rstudio.org). Basic R is essential, but RStudio will make your life a lot easier, I promise!
3. Install the [synapse python client](https://python-docs.synapse.org/build/html/index.html), create a [`.synapseConfig` file](https://python-docs.synapse.org/build/html/Credentials.html) in your home directory.
4. Clone this repository - it has all that you will need to contribute and run this analysis.

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
| Metadata file | [syn25807733](https://www.synapse.org/#!Synapse:1syn25807733) |

The data was pushed from raw files to long-form tables for facile querying and viewing:

| Description | Normalization/filtering| Link |
| --- | --- | --- |
| Global Proteomics | Uncorrected |[syn25808625](https://www.synapse.org/#!Synapse:syn25808625) |
| Global Proteomics | Batch-corrected, no missing batches | [syn25808020](https://www.synapse.org/#!Synapse:syn25808020)|
| Global phosphoproteomics | Uncorrected |[syn25808685](https://www.synapse.org/#!Synapse:syn25808685)|
| Global phosphoprotoemics | Batch-corrected, no missing batches |[syn25808662](https://www.synapse.org/#!Synapse:syn25808662)|
| Global phosphoproteomics | Batch-corrected, at most 10 missing batches |[syn26469873](https://www.synapse.org/#!Synapse:syn26469873)|
| Global phosphoproteomics | Batch-corrected, at most 4 missing batches |[syn26477193](https://www.synapse.org/#!Synapse:syn26477193) |

Before batch correction, we filter to eliminate features which contain too much missing data. Originally, we filtered out any features which had at least one batch with ONLY missing data in that batch, ie, "no missing batches". However, this filtering was too strict for our large phosphoproteomics dataset, so we applied two less conservative filters prior to batch correction, as the "no missing batches" filter left us with relatively few phosphosites. Thankfully, as our global proteomics dataset has very little missing data, the same issue did not arise here.

### Gene mutation, RNASeq, Clinical data

These data are also uploaded to Synapse, and then parsed.

| Description | Link |
| --- | -- |
| Waves 1to4 WES Data | [syn2648827](https://www.synapse.org/#!Synapse:syn26428827/tables/) |
| Waves 1to4 RNASeq Data | [syn26428813](https://www.synapse.org/#!Synapse:syn26428813) |
| Drug response data | [syn25813252](https://www.synapse.org/#!Synapse:syn25830473)|


These tables were pulled from spreadsheets that are currently also stored on Synapse in [this folder](https://www.synapse.org/#!Synapse:syn24171152). Clinical data is still being updated, but is currently stored in an [excel spreadsheet](https://www.synapse.org/#!Synapse:syn25796769).

## Data Analysis

This section aspirationally aims to serve as the outline for the manuscript we are building. This is pretty rough so as things merge together or overlap we can consider restructuring the repository to reflect the latest state of the manuscript.

### Cohort exploration and summarization

The first figure of the manuscript will require visualizing the Beat AML cohort and the data we have. To date this analysis requires, the following, each of which should produce either data for future analysis or figure panels for figure 1. The code should be deposited in the [cohort_summary/](./cohort_summary) directory.

#### Circos plot of data types

We are looking into circos plotting to summarize the data types, this will enable us to see how much data there are for each patient.

#### Non-negative matrix factorization

This will take a multi-omics approach to clustering all samples. This clusters and metagenes will be stored on synapse for future anlaysis.

#### Patient stratification by cluster

Now that we have the patient assignments we can ask if there are survival differences between patients of each cluster, if there are genetic mutation differences, or if there are other clinical properties that vary.

#### Metagene analysis

Last we need to investigate the 'metagenes' that define the clusters and determine if there is functional enrichment, or any phospho networks that are active or depleted.

### Mutational profiling

Figure 2 will compare the genetic mutation data to the other data types. This code will go into the [mutational_analysis](./mutational_analysis/) directory. This analysis focuses on the transcriptomic and proteomic differences between patients with various mutational combinations.

### Immune deconvolution
Figure 3 will focus on the immune infiltration of various tumors using the BayesDeBulk analysis within the Decomprolute framework.

### Functional analysis
Lastly we will spend Figures 4 and 5 investigating the drug profiles that are unique to this dataset.

#### Regression derived signatures
What can we learn from the regression? Are there gene sets of interest?

#### Differential expression
Should we do differential expression as well?

## Paper Figures
This section is reserved for the manuscript figures, as we work on the analysis.
