![Bionano logo](images/Bionano-Logo.png?raw=true)

## Enrichment Project Workflow Development
---
Author: Syukri Shukor

Contact: sshukor@bionano.com or sykrishukor@gmail.com

This repository documents python notebooks and scripts used for Clinical Affairs cohort studies, as well as various applications that require manipulation and analysis of SMAP files. More information on SMAP files can be found in [OGM File Format Specification Sheet, Page 30](https://bionano.com/wp-content/uploads/2023/08/CG-00045-OGM-File-Format-Specification-Sheet.pdf)

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

Here are highlights of functions and descriptions.

Remove # headers from SMAP file, keeping the column headers.
```
get_smap_table(smap_path, output_path)
```


Excludes SV types that are false positive.
```
exclude_types_svs(df, nonTrans=False)
```


Filter SV dataframe according to confidence scores. Edit confidence score cutoffs according to what is recommended by Access in a sample's circos plot.
```
filter_sig_svs(df)
```


'Flattens' inversions with its pairs by matching the SmapID with LinkID.
```
flatten_df_inv(df)
```


Finds genes in bed_df that overlap an SV breakpoints.
```
filter_genes(bed_df, chr, start, end)
```


Standardizes output of `filter_genes` function by sorting genes in alphabetical order.
```
parse_sorted_unique_genes(gene_row)
```


Remaps SV types to its 3 character designation.
```
get_sv_type(type_str)
```


Cluster TWO SV dataframes into one by 'OverlapGenes', 'chr', 'Type', 'Zygosity', 'RefStartPos',and 'RefEndPos', grouping similar SVs into a cluster by specified `posWindow `and `reciprocalSize`. By default `posWindow` is 5kb and `reciprocalSize` is 50%.

```
# main driver function
group_sv_by_gene_dual_cohort(cleaned_df)

# group_sv_by_gene_dual_cohort calls the following function
cluster_refpos_dual(prev_cluster_ID, df_in, posWindow, reciprocalSize)
```


Cluster ONE SV dataframe by 'OverlapGenes', 'chr', 'Type', 'Zygosity', 'RefStartPos',and 'RefEndPos', grouping similar SVs into a cluster by specified `posWindow `and `reciprocalSize`. By default `posWindow` is 5kb and `reciprocalSize` is 50%.

```
# main driver function
group_sv_by_gene_single_cohort(cleaned_df)

# group_sv_by_gene_single_cohort calls the following function
cluster_refpos_single(prev_cluster_ID, df_in, posWindow, reciprocalSize)
```


Cluster ONE SV dataframe by 'chr', 'Type', 'Zygosity', 'RefStartPos',and 'RefEndPos', grouping similar SVs into a cluster by specified `posWindow `and `reciprocalSize`. By default `posWindow` is 5kb and `reciprocalSize` is 50%.

```
# main driver function
group_sv_single_cohort_nogenes(cleaned_df)

# cluster_sv_single_cohort_nogenes also calls the following function
cluster_refpos_single_nogenes(prev_cluster_ID, df_in, posWindow, reciprocalSize)
```

## EXAMPLE USAGE

See `~/jupyter_notebooks` for .ipynb that uses variations of the scripts. Select `smap_analysis` as the kernel before running any code in the notebooks.

Import the util py script by declaring:
```
from smap_analysis_util import *
```

For example, `partial_dup_vs_ctrldb_cluster_unique_id_overlap_exons.ipynb` uses most of the functions above including.