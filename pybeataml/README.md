# TODO code/style

- Write load_beatAML_data.R using the data class.
- Create pipeline for enrichment from gene sets.
- Connect enriched terms to AUC predictions
- Write documentation to allow others to setup and test.
- ~~Create data class to load all in a pretty format.~~
- ~~Create/upload sample model code in python.~~

# TODO analysis

- Create merged/hybrid enrichment analysis of various data types
- Generate networks from enriched terms using MAGINE

# Install miniconda and create env

Download and install miniconda (if you'd like to use the one liner, if not you can install the libraries in
environment.yml any which way you'd like, python major flaw = too many options).

`conda env create -f environment.yml`

`conda activate beatAML_env_37`

Strong likelihood one needs to clone the latest version of MAGINE locally and add to $PYTHONPATH. Will release new
MAGINE version when I have more time.

# Data class

There is a single class interface for the BeatAML datasets. This can be found in `load_data.py`. A few things to note.
All data has a gene_symbol and label column. The gene_symbol is just that, where the label will have additional
information depending on the data source. For instance, rna_seq will have the gene symbol plus `_rna` appended, so when
one trains a model all genes can be separated with names rather than `BAX.1`. For the WES data, there are two
level : `wes` which uses gene summarized binary or `wes_protein_level`
which keeps the amino acid number and letter change (`BAX-A152G`). Could be useful to not drop information for later
analysis.

``` python
from pybeataml.load_data import AMLData
data = AMLData()
# get trainable dataset in patient x feature format
source = ['proteomics', 'rnaseq'] # subset of data types
drug_name = 'Sorafenib' 
subset = data.get_trainable_data(source, drug_name)

# get specific data sets from class (simple convention for attributes)
mutations = data.wes
proteomics = data.proteomics
phospho = data.phospho
rna = data.rna

#subset data by datatype
# Long table format
rna_protein = data.subset_flat(['rnaseq', 'proteomics'])
# pivoted to be sample x label
rna_protein = data.subset(['rnaseq', 'proteomics'])
```

# Examples of using python calls in R

- `example_synapse_via_python.R`
- `example_python_dataclass_in_R.R`

```R
library(reticulate)
reticulate::use_condaenv("beatAML_env_37", required = TRUE)
# might need to adjust path
reticulate::source_python("load_data.py")

# load data class
data <- AMLData()

# can access functions with $
mat <- data$all_data

# same as python examples above, but use $
phospho <- data$phospho
rna <- data$rna

# function calls are the same as above as well
rna_protein <- data.subset(['rnaseq', 'proteomics'])

```

# Spin up Jupyter notebook

`jupyter notebook`

# Summary notebooks

- `run_models.ipynb`
  - Development notebook to create models. Most code starts here then is refactored to .py files
- `meta_gene_enrichment.ipynb`
  - Enrichment analysis on meta genes from NMF
- `plot_model_summaries.ipynb`
  - Visualize outputs from models/run_all_kfold_ultimate.py
# Regression models

A single script to run through Lasso, ElasticNet, Gradient Boosted Trees, and Support Vector Regression can be found
in `models/run_all_kfold_ultimate.py`. It runs through all drugs, all data types individually, and in combinations.
Typically, run time takes > 1 hour, but there are a few threading flags that can be improved. In the same folder there
are a few variations of similar scripts. `run_all.py` is a simpler script, but doesn't output features selected from the
modeling, however is quite faster and easier to follow. The other two scripts are both test beds for scanning
parameters. Rather than run k-fold over all combinations of drugs, data types, and models, I thought it might be cleaner
to follow up on optimizing models after they are identified.

