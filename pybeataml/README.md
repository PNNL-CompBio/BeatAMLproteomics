# In progress

Still need to do a lot. Just setting up templates for the goals.

# TODO

- Create data class to load all in a pretty format.
- Create pipeline for enrichment from gene sets.
- Connect enriched terms to AUC predictions
- Write up documentation to allow others to setup and test.
- Create/upload sample model code in python.
- Compare R and python model predictions.

# Install miniconda

`conda env create -n environment.yml`

`conda activate beatAML_env_37`

Strong likelikehood one needs to clone the latest version of MAGINE locally and add to $PYTHONPATH. Will release new
MAGINE version when I have more time.

# Spin up Jupyter notebook

`jupyter notebook`

# Summary notebooks

- run_models.ipynb
    - Development notebook to create models. Most code starts here then is refactored to .py files
- meta_gene_enrichment.ipynb
    - Enrichment analysis on meta genes from NMF

# Examples of using python calls in R

- example_synapse_via_python.R
- example_python_dataclass_in_R.R