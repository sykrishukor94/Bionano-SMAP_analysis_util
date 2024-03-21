![Bionano logo](images/Bionano-Logo.png?raw=true)

## Enrichment Project Workflow Development
---
Author: Syukri Shukor

Contact: sshukor@bionano.com, achitsazan@bionano.com or apang@bionano.com

This repository documents python notebooks and scripts used for Clinical Affairs cohort study. 

## SETUP
---
Current implementation of this workflow utilizes the following packages.
```
conda (4.12.0)
mamba (0.15.3)
Python (3.10.6)
Solve (3.7.2)
```

**Conda installation**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#Follow-installation prompts
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

followed by installing the packages listed in `smap_analysis.yml`

```
conda env create -f smap_analysis.yml
```
This will install jupyterlab, pandas, and seaborn into an isolated software environment, that has to be activated with

```
$ conda activate smap_analysis
```

## SMAP ANALYSIS FUNCTIONS



## RUNNING PYTHON SCRIPTS
---

```
screen
```

```
$(base) conda activate smap_analysis
$(smap_analysis)
```
Execute example script
```
python .\get_percent_presence_ctrl.py 
```
Detach screen
```
ctrl-a ctrl-d
```

For the jupyter notebooks, open .ipynb in VS code. Select `smap_analysis` as the kernel before running any code in the notebooks.