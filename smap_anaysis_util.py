# script compiling enrichment project helper functions

# load packages
from itertools import groupby
# from matplotlib_venn import venn2, venn2_circles
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy

from scipy import stats
from scipy.stats import chi2
from scipy.stats import chisquare
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
from scipy.stats import f_oneway
# from scipy.stats.contingency import odds_ratio

import copy
import re 
import os
import pathlib


# Standardized headers for both bed files
# headers: chrom	chromStart	chromEnd	Gene	Index	strand	chromStart2	chromEnd	RGB
bed_headers = ["chr","RefStartPos","RefEndPos","Gene","Index","Strand","RefStartPos2","RefEndPos2","RGB"]

# Load hg38_gaps.bed
gaps_df = pd.read_csv(r"/home/users6/sshukor/enrichment_project/input/ref_genes/hg38_gaps.bed", dtype=object,  index_col=False, sep="\t", names=bed_headers)
gaps_df = gaps_df.astype({'RefStartPos':'float64', 'RefEndPos':'float64', 'chr':'int'})



# Helper function to get smap table without preceding "#..."" lines
def get_smap_table(smap_path, output_path):

    output = open(f"{output_path}", "w")

    with open(smap_path) as file:
        for line in file:

            if line[0:2] == "#h":
                output.write(line[3:])
            elif line[0:1] == "#":
                pass
            else:
                output.write(line)

    output.flush()

    # "reset" fd to the beginning of the file
    output.seek(0)
    # print("File contents:\n", output.read())

    print(f"output file path: {output_path}.smap")
    output.close()


# filter_rare_svs(df, percent_cutoff)​
# Helper function to filter for SVs equal or below a certain threshold (inadvertently excludes inversion partials)
def filter_rare_svs(df, percent_cutoff):
    df = df.astype({"Present_in_%_of_BNG_control_samples":'float64', "Present_in_%_of_BNG_control_samples_with_the_same_enzyme":'float64'})

    return df.loc[((df["Present_in_%_of_BNG_control_samples_with_the_same_enzyme"] <= percent_cutoff))
            & ((df["Present_in_%_of_BNG_control_samples_with_the_same_enzyme"] != -1))]


# exclude_types_svs(df)​
def exclude_types_svs(df, nonTrans=False):
    if nonTrans == True:
        excluded_types = ["deletion_nbase", "insertion_nbase", "gain_masked", "loss_masked", "inversion_partial", "trans_interchr_common", "trans_intrachr_common", 'translocation_interchr', 'trans_intrachr_segdupe', 'trans_interchr_segdupe', 'translocation_intrachr', 'insertion_tiny', 'deletion_tiny']
    else:
        excluded_types = ["deletion_nbase", "insertion_nbase", "inversion_nbase" "gain_masked", "loss_masked", "inversion_partial", "trans_interchr_common", "trans_intrachr_common", 'trans_intrachr_segdupe', 'trans_interchr_segdupe', 'insertion_tiny', 'deletion_tiny']
    return df.loc[(~df["Type"].isin(excluded_types)) | df['Type'].str.contains('')]


# ['deletion', 'insertion', 'insertion_nbase',
#        'trans_interchr_common', 'trans_intrachr_common',
#        'translocation_interchr', 'duplication', 'inversion',
#        'trans_intrachr_segdupe', 'duplication_inverted', 'deletion_nbase',
#        'inversion_paired', 'trans_interchr_segdupe', 'duplication_split',
#        'translocation_intrachr']

# Helper function to filter by confidence intervals

# 'trans', 'duplication_inverted', 'duplication','duplication_split', 'inv', 'ins', 'del'
def filter_sig_svs(df):
    dict_conf = {
        "ins": 0,
        "del": 0,
        "dup": -1,
        "inv": float(0.7),
        "trans": float(0.05) # intrachr and interchr
    }

    def get_sv_conf_thres(sv):
        sv = str(sv)
        if "ins" in sv:
            return dict_conf["ins"]
        elif "del" in sv:
            return dict_conf["del"]
        elif "duplication" in sv:
            return dict_conf["dup"]
        elif "inv" in sv:
            return dict_conf["inv"]
        elif "trans" in sv:
            return dict_conf["trans"]
        else:
            print(f"unknown sv type: {sv}")

    df = df.astype({'Confidence':'float64'})
    pass_conf = df[["Confidence", "Type"]].apply(lambda x: True if x[0] >= float(get_sv_conf_thres(x[1])) else False, axis=1)

    return df.loc[pass_conf == True]


# flatten_df_inv(df)​
def flatten_df_inv(df):
    # print(df.columns)
    print("Total:", len(df))

    # subset inversions and inversion_partials
    # inv = ['inversion', 'inversion_partial', 'inversion_repeat']
    inv_df = (df.loc[df['Type'].str.contains('inversion')]).copy()

    print("inversions, inversion nbase, and inversion & inversion_partials:", len(inv_df))

    # initialize empty dataframe
    out_df = pd.DataFrame()
    # proceed if there are inversions
    if len(inv_df) > 0:
        # remove inversions from main df
        df = df.loc[~df['Type'].str.contains('inversion')]
        # print(df.head())

        # initialize column of T/F, where SVs == T have been analyzed for flattening, and F haven't
        inv_df.loc[:, 'flattened'] =  False

        # loops through inversions and flattens inversions where LinkID matches SmapID
        while False in pd.unique(inv_df['flattened']):
            
            df_copy = copy.deepcopy(inv_df.loc[inv_df['flattened'] == False])
            # print(len(df_copy))

            # get first row item
            qry_ID = df_copy['Sample_ID'].iat[0]
            qry_chr = df_copy['chr'].iat[0]
            qry_smapID = df_copy['SmapID'].iat[0]
            qry_linkID = df_copy['LinkID'].iat[0]

            # find pair for the row if exists
            df_copy = df_copy.loc[(df_copy['Sample_ID'] == qry_ID)
                                & (df_copy['chr'] == qry_chr)
                                & ((df_copy['SmapID'] == qry_linkID)
                                | (df_copy['SmapID'] == qry_smapID))]

            # update original inv_df column as flattened
            inv_df.loc[(inv_df['Sample_ID'] == qry_ID) 
                     & (inv_df['chr'] == qry_chr) 
                     & ((inv_df['SmapID'] == qry_linkID) 
                     | (inv_df['SmapID'] == qry_smapID)), ['flattened']] = True

            # flatten if there is a pair, report only row if not
            if len(df_copy) > 1:

                # choose inversion in the pair, combine both
                if "inversion" in pd.unique(df_copy['Type']):
                    row = df_copy.loc[(df_copy['RefEndPos'] != -1.0)].reset_index() # choose inversions
                else:
                    row = df_copy.iloc[[0]].reset_index()
                
                # merge pairs, choosing the lowest and largest
                ref = np.concatenate((df_copy['RefStartPos'].to_numpy(), df_copy['RefEndPos'].to_numpy()), axis=None)
                ref = np.delete(ref, np.where(ref == -1))
                # print(ref)

                row.at[0, 'RefStartPos'] = np.amin(ref, axis=None)
                row.at[0, 'RefEndPos'] = np.amax(ref, axis=None)
                row.at[0, 'SVsize'] = row['RefEndPos'] - row['RefStartPos']

            else:
                row = df_copy.iloc[[0]].reset_index()
            # append to out_df        
            # print(row)
            out_df = pd.concat([out_df, row], axis=0, ignore_index=True)

    # append out_df to df

    # print(df.columns)
    return pd.concat([out_df[df.columns], df], axis=0, ignore_index=True)


# filter_genes(bed_df, chr, start, end)​
# Helper function to remove SVs overlapping gaps bed file
def filter_genes(bed_df, chr, start, end):
    
    match_chr = bed_df.loc[(bed_df["chr"] == chr)]

    # display(match_chr)

    # transform matched chr rows into numpy (faster processing than directly panda slicing)
    bed_genes = np.array(match_chr["Gene"])
    bed_start = np.array(match_chr["RefStartPos"])
    bed_end = np.array(match_chr["RefEndPos"])

    overlap_idx = np.where((((start <= bed_start) & (bed_start <= end)) 
                            | ((bed_start <= start) & (start <= bed_end))))
    
    overlap_genes = bed_genes[overlap_idx]
    
    if len(overlap_genes) > 0:
        return ";".join(list(np.unique(overlap_genes)))
    else:
        return "-"
    

def parse_sorted_unique_genes(gene_row):
    # return gene_col.str.split(';').map(lambda x: ';'.join(list(set(sorted(x)))))
    gene_row = list(set(gene_row.split(';')))
    
    return ';'.join(sorted(gene_row))


# Helper function to filter genes
# Helper function to remove SVs overlapping gaps bed file
def filter_nogaps(type, chr, start, end):

    # n_gaps check only on insertion and deletions
    if (("ins" in type) | ("del" in type)): 
        match_chr = gaps_df.loc[(gaps_df["chr"] == chr)]
        # check each SV related to gaps[chr"] == chr for overlaps, False if none

        # transform matched chr rows into numpy (faster processing than directly panda slicing)
        # bed_genes = np.array(match_chr["Gene"])
        bed_start = np.array(match_chr["RefStartPos"])
        bed_end = np.array(match_chr["RefEndPos"])

        for bedStart, bedEnd in zip(bed_start, bed_end):
            if (((start <= bedStart) & (bedStart <= end)) 
            | ((bedStart <= start) & (start <= bedEnd))):
                return True

    return False


# Helper function to re-map SV types
def get_sv_type(type_str):
    # type_str = type_str.strip()
    if "ins" in type_str:
        return "ins"
    elif "del" in type_str:
        return "del"
    elif "duplication" in type_str:
         # ignore duplications to stay in concordance with ctrl_df['Type]
        return type_str
    elif "inv" in type_str:
        return "inv"
    elif "trans" in type_str:
        return "trans"
    else:
        return(type_str)


# Helper function to cluster
# 2 Helper functions to cluster SVs. Clustering removes redundant SVs, keeping only `unique` SVs
pd.options.mode.chained_assignment = None
# Helper function to cluster gene coordinates, given a group of SV(s) grouped by gene, chr, zygo, and type
def cluster_refpos_dual(prev_cluster_ID, df_in, posWindow, reciprocalSize):
    
    # Get sv info from first df item
    sv_info_headers = ["OverlapGenes", "chr", "Type", 'Zygosity']
    sv_info = df_in[sv_info_headers].head(1)

    qry_type = df_in["Type"].iat[0]

    sv_df = pd.DataFrame()
    df_out = pd.DataFrame()
    while False in pd.unique(df_in['clustered']):
        
        # subset non-clustered SVs and get first item
        # df_in = df_in.loc[df_in['clustered'] == False, ['RefStartPos', 'RefEndPos', 'SVsize', 'case_ID', 'ctrl_ID', 'num_overlap_DGV_calls', 'clustered']]
        df_in = df_in.loc[df_in['clustered'] == False]
        
        # get gene coordinates and size of first item
        qry_start = float(df_in['RefStartPos'].iat[0])
        qry_end = float(df_in['RefEndPos'].iat[0])
        qry_size = float(df_in['SVsize'].iat[0])

        # find overlaps with first item (special treatnment for inversions since no inv size in ctrldb)
        if qry_type in ['ins','insertion', 'del' ,'deletion', 'dup','duplication_paired' ,'duplication_inverted', 'duplication', 'duplication_split']:
        
            cluster_df = df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))) \
                                    & (((qry_size/df_in['SVsize'])*100 >= reciprocalSize) \
                                        & ((df_in['SVsize']/qry_size)*100 >= reciprocalSize))]
            
            # mark as clustered in original input_df
            df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))) \
                                    & (((qry_size/df_in['SVsize'])*100 >= reciprocalSize) \
                                        & ((df_in['SVsize']/qry_size)*100 >= reciprocalSize)), 'clustered'] = True

        else: # else if its inv (or trans)

            cluster_df = df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos'])))]
            
            df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))), 'clustered'] = True

        # get cluster info
        # get SV key information (OverlapGenes, Chr, Type, Zygo)
        row = sv_info
        
        # set cluster and row ID number
        prev_cluster_ID += 1
        cluster_df.loc[:, 'cluster_ID'] = prev_cluster_ID
        row['cluster_ID'] = prev_cluster_ID

        # get list of gene coordinates 
        row['RefStartPos'] = ", ".join(cluster_df['RefStartPos'].astype(str))
        row['RefEndPos'] = ", ".join(cluster_df['RefEndPos'].astype(str))

        # Get info on n samples that share an SV, ignoring repeats
        len_case = len(cluster_df['case_ID'].dropna().unique())
        len_ctrl = len(cluster_df['ctrl_ID'].dropna().unique())

        row['num_SVs'] = len(cluster_df)

        row['total_counts'] = len_case + len_ctrl
        row['case_counts'] = len_case    # get # of unique samples
        row['ctrl_counts'] = len_ctrl    # get # of unique samples

        # process sample IDs (include repeats) into a list
        case_id = cluster_df['case_ID'].to_numpy()
        ctrl_id = cluster_df['ctrl_ID'].to_numpy()
        
        ids = []
        for i in range(len(case_id)):
            if pd.isna(case_id[i]):
                ids.append(ctrl_id[i])
            else:
                ids.append(case_id[i])

        row['Sample_ID'] =  ", ".join(ids)

        # populate lists of other sample-specific metrics for table QC

        row['SVsize'] = ", ".join(cluster_df['SVsize'].astype(str))
        row['sv_means'] = cluster_df['SVsize'].mean()
        row['sv_std'] = cluster_df['SVsize'].std()
        row['case_oDgv'] = ", ".join(list(pd.unique(cluster_df['num_overlap_DGV_calls'].dropna())))
        # print(row)

        # remove redundant SVs
        cluster_df['merged ID'] = cluster_df.apply(lambda x: x['ctrl_ID'] if pd.isna(x['case_ID']) else x['case_ID'], axis=1)
        cluster_df = cluster_df.drop_duplicates(subset=['merged ID'])


        sv_df = pd.concat([cluster_df, sv_df], axis=0, ignore_index=True)
        df_out = pd.concat([df_out, row], axis=0, ignore_index=True)
        
    return prev_cluster_ID, sv_df, df_out

def group_sv_by_gene_dual_cohort(cleaned_df):
    # initialize column of T/F, where SVs == T have been clustered, and F haven't.
    cleaned_df.loc[:, 'clustered'] = False    
    cleaned_df.loc[:, 'cluster_ID'] = -1

    # Sort first by refstartstop from smallest to largest
    cleaned_df = cleaned_df.sort_values(["OverlapGenes", "chr", "Type", 'Zygosity', 'RefStartPos','RefEndPos'])

    # Group compiled rare SVs by criteria, then list samples sharing each rare SV.
    # due to memory issues, we put df_g in a dictionary, using identified rare SVs as key, then iterate by key to values containing df
    groups = dict(list(cleaned_df.groupby(by=["OverlapGenes", "chr", "Type", 'Zygosity'])))

    df_in = pd.DataFrame()
    df_out = pd.DataFrame()

    total = len(groups.values())
    count = 0
    prev_cluster_ID = 0

    # main loop to process every SV group into smaller clusters, appends clustered dataframe into out_df
    for df in groups.values():
        prev_cluster_ID, cluster_df_in, cluster_df_out = cluster_refpos_dual(prev_cluster_ID, df, 5000.0, 50)

        # returned cluster SVs are appended to output df.
        df_in = pd.concat([df_in, cluster_df_in], axis=0, ignore_index=True)
        df_out = pd.concat([df_out, cluster_df_out], axis=0, ignore_index=True)

        count += 1
        if count % (total/10) == 0:
            print(f"{(count/total)*100}%/ 100% done")

    return df_in, df_out



# Single cohort clustering
# 2 Helper functions to cluster SVs. Clustering removes redundant SVs, keeping only `unique` SVs
pd.options.mode.chained_assignment = None
# Helper function to cluster gene coordinates, given a group of SV(s) grouped by gene, chr, zygo, and type
def cluster_refpos_single(prev_cluster_ID, df_in, posWindow, reciprocalSize):
    
    # Get sv info from first df item
    sv_info_headers = ["OverlapGenes", "chr", "Type", 'Zygosity']
    sv_info = df_in[sv_info_headers].head(1)

    qry_type = df_in["Type"].iat[0]

    sv_df = pd.DataFrame()
    row_df = pd.DataFrame()
    while False in pd.unique(df_in['clustered']):
        
        # subset non-clustered SVs and get first item
        # df_in = df_in.loc[df_in['clustered'] == False, ['RefStartPos', 'RefEndPos', 'SVsize', 'case_ID', 'ctrl_ID', 'num_overlap_DGV_calls', 'clustered']]
        df_in = df_in.loc[df_in['clustered'] == False]
        
        # get gene coordinates and size of first item
        qry_start = float(df_in['RefStartPos'].iat[0])
        qry_end = float(df_in['RefEndPos'].iat[0])
        qry_size = float(df_in['SVsize'].iat[0])

        # find overlaps with first item (special treatnment for inversions since no inv size in ctrldb)
        if qry_type in ['ins','insertion', 'del' ,'deletion', 'dup','duplication_paired' ,'duplication_inverted', 'duplication', 'duplication_split']:
        
            cluster_df = df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))) \
                                    & (((qry_size/df_in['SVsize'])*100 >= reciprocalSize) \
                                        & ((df_in['SVsize']/qry_size)*100 >= reciprocalSize))]
    
            
            # mark as clustered in original input_df
            df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))) \
                                    & (((qry_size/df_in['SVsize'])*100 >= reciprocalSize) \
                                        & ((df_in['SVsize']/qry_size)*100 >= reciprocalSize)), ['clustered']] = True

        else: # else if its inv (or trans)

            cluster_df = df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos'])))]
            
            df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))), ['clustered']] = True


        # get cluster info
        # get SV key information (OverlapGenes, Chr, Type, Zygo)
        row = sv_info
        
        # set cluster and row ID number
        prev_cluster_ID += 1
        cluster_df['cluster_ID'] = prev_cluster_ID
        row['cluster_ID'] = prev_cluster_ID

        # get list of gene coordinates 
        row['RefStartPos'] = cluster_df['RefStartPos'].min()
        row['RefEndPos'] = cluster_df['RefEndPos'].max()
        row['listRefStartPos'] = ", ".join(cluster_df['RefStartPos'].astype(str))
        row['listRefEndPos'] = ", ".join(cluster_df['RefEndPos'].astype(str))

        # Get info on n samples that share an SV, including and excluding repeats
        row['num_SVs'] = len(cluster_df)
        row['num_unique_samples'] = len(cluster_df['Sample_ID'].dropna().unique())

        row['Sample_ID'] =  ", ".join(cluster_df['Sample_ID'])

        # populate lists of other sample-specific metrics for table QC

        row['SVsize'] = ", ".join(cluster_df['SVsize'].astype(str))
        row['sv_means'] = cluster_df['SVsize'].mean()
        row['sv_std'] = cluster_df['SVsize'].std()
        # row['num_overlap_DGV_calls'] = ", ".join(list(pd.unique(cluster_df['num_overlap_DGV_calls'].dropna())))
        # print(row)

        # display(cluster_df)
        # Remove duplicate sample_ID and keep only 1 SV in a cluster per sample
        # cluster_df = cluster_df.drop_duplicates(subset=['Sample_ID'])
        # display(cluster_df)
        
        sv_df = pd.concat([cluster_df, sv_df], axis=0, ignore_index=True)
        row_df = pd.concat([row_df, row], axis=0, ignore_index=True)
        
    return prev_cluster_ID, sv_df, row_df


def group_sv_by_gene_single_cohort(cleaned_df):
    # initialize column of T/F, where SVs == T have been clustered, and F haven't.
    cleaned_df.loc[:, 'clustered'] = False    
    cleaned_df.loc[:, 'cluster_ID'] = -1

    # Sort first by refstartstop from smallest to largest
    cleaned_df = cleaned_df.sort_values(["OverlapGenes", "chr", "Type", 'Zygosity', 'RefStartPos','RefEndPos'])

    # Group compiled rare SVs by criteria, then list samples sharing each rare SV.
    # due to memory issues, we put df_g in a dictionary, using identified rare SVs as key, then iterate by key to values containing df
    groups = dict(list(cleaned_df.groupby(by=["OverlapGenes", "chr", "Type", 'Zygosity'])))

    df_in = pd.DataFrame()
    df_out = pd.DataFrame()

    total = len(groups.values())
    prev_cluster_ID = 0

    # main loop to process every SV group into smaller clusters, appends clustered dataframe into out_df
    for df in groups.values():
        # inputs certain columns for clustering
        prev_cluster_ID, cluster_df_in, cluster_df_out = cluster_refpos_single(prev_cluster_ID, df, 5000.0, 50)

        # returned cluster SVs are appended to output df.
        df_in = pd.concat([df_in, cluster_df_in], axis=0, ignore_index=True)
        df_out = pd.concat([df_out, cluster_df_out], axis=0, ignore_index=True)

        print(prev_cluster_ID)
        
        if prev_cluster_ID % (total/10) == 0:
            print(f"{(prev_cluster_ID/total)*100}%/ 100% done")

    return df_in, df_out


# Single cohort clustering
# 2 Helper functions to cluster SVs. Clustering removes redundant SVs, keeping only `unique` SVs
pd.options.mode.chained_assignment = None
# Helper function to cluster gene coordinates, given a group of SV(s) grouped by gene, chr, zygo, and type
def cluster_refpos_single_nogenes(prev_cluster_ID, df_in, posWindow, reciprocalSize):
    
    # Get sv info from first df item
    sv_info_headers = ["OverlapGenes","Zygosity","chr", "Type"]
    sv_info = df_in[sv_info_headers].head(1)
    # display(df_in.head(1))

    qry_type = df_in["Type"].iat[0]

    sv_df = pd.DataFrame()
    row_df = pd.DataFrame()
    while False in pd.unique(df_in['clustered']):
        
        # subset non-clustered SVs and get first item
        # df_in = df_in.loc[df_in['clustered'] == False, ['RefStartPos', 'RefEndPos', 'SVsize', 'case_ID', 'ctrl_ID', 'num_overlap_DGV_calls', 'clustered']]
        df_in = df_in.loc[df_in['clustered'] == False]
        
        # get gene coordinates and size of first item
        qry_start = float(df_in['RefStartPos'].iat[0]) + posWindow
        qry_end = float(df_in['RefEndPos'].iat[0]) + posWindow
        qry_size = float(df_in['SVsize'].iat[0])

        # find overlaps with first item (special treatnment for inversions since no inv size in ctrldb)
        if qry_type in ['ins','insertion', 'del' ,'deletion', 'dup','duplication_paired' ,'duplication_inverted', 'duplication', 'duplication_split']:
        
            cluster_df = df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))) \
                                    & (((qry_size/df_in['SVsize'])*100 >= reciprocalSize) \
                                        & ((df_in['SVsize']/qry_size)*100 >= reciprocalSize))]
    
            
            # mark as clustered in original input_df
            df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))) \
                                    & (((qry_size/df_in['SVsize'])*100 >= reciprocalSize) \
                                        & ((df_in['SVsize']/qry_size)*100 >= reciprocalSize)), 'clustered'] = True

        else: # else if its inv (or trans)

            cluster_df = df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos'])))]
            
            df_in.loc[(((qry_start - posWindow <= df_in['RefStartPos']) & (df_in['RefStartPos'] <= qry_end + posWindow)) \
                                    | ((df_in['RefStartPos'] <= qry_start - posWindow) & (qry_start - posWindow <= df_in['RefEndPos']))), 'clustered'] = True


        # get cluster info
        # get SV key information (OverlapGenes, Chr, Type, Zygo)
        row = sv_info
        
        # set cluster and row ID number
        prev_cluster_ID += 1
        cluster_df['cluster_ID'] = prev_cluster_ID
        row['cluster_ID'] = prev_cluster_ID

        # get list of gene coordinates 
        row['RefStartPos'] = cluster_df['RefStartPos'].min()
        row['RefEndPos'] = cluster_df['RefEndPos'].max()
        row['listRefStartPos'] = ", ".join(cluster_df['RefStartPos'].astype(str))
        row['listRefEndPos'] = ", ".join(cluster_df['RefEndPos'].astype(str))

        # Get info on n samples that share an SV, including and excluding repeats
        row['num_SVs'] = len(cluster_df)
        row['num_unique_samples'] = len(cluster_df['Sample_ID'].dropna().unique())

        row['Sample_ID'] =  ", ".join(cluster_df['Sample_ID'])

        # populate lists of other sample-specific metrics for table QC

        row['SVsize'] = ", ".join(cluster_df['SVsize'].astype(str))
        row['sv_means'] = cluster_df['SVsize'].mean()
        row['sv_std'] = cluster_df['SVsize'].std()
        # row['num_overlap_DGV_calls'] = ", ".join(list(pd.unique(cluster_df['num_overlap_DGV_calls'].dropna())))
        # print(row)

        # display(cluster_df)
        # Remove duplicate sample_ID and keep only 1 SV in a cluster per sample
        # cluster_df = cluster_df.drop_duplicates(subset=['Sample_ID'])
        cluster_df['clustered'] = True
        
        sv_df = pd.concat([cluster_df, sv_df], axis=0, ignore_index=True)
        row_df = pd.concat([row_df, row], axis=0, ignore_index=True)
        
    return prev_cluster_ID, sv_df, row_df


def group_sv_single_cohort_nogenes(cleaned_df):
    # initialize column of T/F, where SVs == T have been clustered, and F haven't.
    cleaned_df.loc[:, 'clustered'] = False    
    cleaned_df.loc[:, 'cluster_ID'] = -1

    # Sort first by refstartstop from smallest to largest
    cleaned_df = cleaned_df.sort_values(["chr", "Type", 'RefStartPos','RefEndPos', 'SVsize'])

    # Group compiled rare SVs by criteria, then list samples sharing each rare SV.
    # due to memory issues, we put df_g in a dictionary, using identified rare SVs as key, then iterate by key to values containing df
    groups = dict(list(cleaned_df.groupby(by=["chr", "Type"])))

    df_in = pd.DataFrame()
    df_out = pd.DataFrame()

    total = len(groups.values())
    prev_cluster_ID = 0

    # main loop to process every SV group into smaller clusters, appends clustered dataframe into out_df
    for df in groups.values():
        # inputs certain columns for clustering
        prev_cluster_ID, cluster_df_in, cluster_df_out = cluster_refpos_single_nogenes(prev_cluster_ID, df, 5000.0, 50)

        # returned cluster SVs are appended to output df.
        df_in = pd.concat([df_in, cluster_df_in], axis=0, ignore_index=True)
        df_out = pd.concat([df_out, cluster_df_out], axis=0, ignore_index=True)

        print(prev_cluster_ID)

    return df_in, df_out



