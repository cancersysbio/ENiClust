#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on Sep 14, 2022

@author: prisnir, lisemangiante
'''

import os
import pandas as pd
import numpy as np
from sndhdr import test_8svx
from pickle import FALSE
from collections import defaultdict

from eniclust.utils import get_data_dir, create_dir
from eniclust import constants

import argparse

from sklearn import datasets, neighbors, linear_model, ensemble
from sklearn.model_selection import cross_val_score, StratifiedKFold, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score, classification_report
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from itertools import chain

import warnings
warnings.filterwarnings('ignore')

class PreProc:
    '''
    classdocs
    '''

    def __init__(self, gene_copy_number_file, seg_file, facets_qc_file, clin_meta_file, iC10_results_file,
                    cytoband_map = os.path.join(get_data_dir(), 'general/exome_gene_locations_arminfo.txt'), 
                    keep_clusters = ['ic1', 'ic10', 'ic2', 'ic4', 'ic5', 'ic6', 'ic8', 'ic9']):
        '''
        Constructor
        '''
        self.gene_copy_number_file = gene_copy_number_file
        self.seg_file = seg_file
        self.facets_qc_file= facets_qc_file
        self.clin_meta_file = clin_meta_file

        if iC10_results_file is not None:
            self.iC10_results_file = iC10_results_file
        
        self.cytoband_map = cytoband_map
        self.keep_clusters = keep_clusters
        
        self.selected_genes = pd.read_csv(constants.training_gene_sets['WES'], sep = '\t')
        print("Is it the right version of the selected genes list? ", ('HIST2H2BE' in self.selected_genes['genes'].tolist()) & ('LY6L' not in self.selected_genes['genes'].tolist()) & (len(self.selected_genes['genes'].tolist()) == 611))
        
        self.genes_trans = self.selected_genes['genes'].tolist()
        self.arms_trans = ['11q13', '8q24', '8p11', '17q24']
        self.other_features_trans = ['genome_doubled', 'ploidy','hrd_loh', 'fraction_cna', 'avg_breakpoints_chr', 'ER', 'HER2', 'WGS', 'WES']

### For 01. gene CN level features
    def _format_gene_CN(self): #to use if we choose to use it instead of the R script I made for caris
        '''Internal method:
        This method reads the gene CN level
        file from Facets. Converts into a sample x gene format '''
        gene_CN_data = pd.read_csv(self.gene_copy_number_file, sep = "\t")
        gene_CN_data_rel = gene_CN_data[["sample", "gene", "mcn"]]
        gene_CN_data_rel.index = gene_CN_data_rel["gene"]
        sample_name = gene_CN_data["sample"].unique()
        gene_CN_data_rel_dict = gene_CN_data_rel[["mcn"]].to_dict()
        gene_CN_data_rel_dic_reformatted = pd.DataFrame.from_dict(gene_CN_data_rel_dict, orient="index")
        gene_CN_data_rel_dic_reformatted.index = sample_name
        return gene_CN_data_rel_dic_reformatted

    def gene_level_CN(self):
        gene_level_CN = pd.read_csv(self.gene_copy_number_file, sep = '\t')
        gene_level_CN = gene_level_CN.set_index('Sample', drop = True)
        return gene_level_CN

    def get_gene_alias(self, df_CN, genes_ref_file):
        '''
        genes_ref_file (Adapted iC10 genes file):
            Probe_ID
            Gene_symbol
            Ensembl_ID
            Cytoband
            Genomic_location_hg18
            chromosome_name_hg18
            start_position_hg18
            end_position_hg18
            Synonyms_0
            Gene.Chosen
            Genomic_location_hg19
            chromosome_name_hg19
            start_position_hg19
            end_position_hg19
            Genomic_location_hg38
            chromosome_name_hg38
            start_position_hg38
            end_position_hg38
            strand
            Alias
            gencode_seqnames_hg19
            gencode_seqnames_hg38
            gencode_start_hg19
            gencode_start_hg38
            gencode_end_hg19
            gencode_end_hg38
            Cytoband_hg19
            Cytoband_hg38
            Source
        '''
        ref_genedef = pd.read_csv(genes_ref_file, sep = "\t")
        gene_list = ref_genedef['Gene_symbol'].unique().tolist()
        print("Right version of the alias annotation file?", len(gene_list) == 611)
        
        wes_genes = df_CN.columns.tolist()
        missing_genes = sorted([i for i in gene_list if i not in df_CN.columns.tolist()])
        print("Number of ENiClust genes with no CN values found in FACETS with their names:", len(missing_genes),"/",len(gene_list))
        
        df_missing_genes = ref_genedef[ref_genedef['Gene_symbol'].isin(missing_genes)][['Gene_symbol','Alias']].sort_values('Gene_symbol')
        df_missing_genes_no_alias = df_missing_genes[(pd.isnull(df_missing_genes['Alias']))]['Gene_symbol'].tolist()
        if len(df_missing_genes_no_alias) != 0:
            print('No alias found for:', df_missing_genes_no_alias)
        
        df_missing_genes = df_missing_genes[~(pd.isnull(df_missing_genes['Alias']))]
        dict_alias_temp = dict(zip(df_missing_genes['Gene_symbol'], df_missing_genes['Alias']))
        
        defdict_alias = defaultdict(list)
        dict_alias = {}
        missing_value = []
        for k, v in dict_alias_temp.items():
            aliases = v.split(";")
            aliases = [i.strip() for i in aliases]
            defdict_alias[k] = aliases
            aliases_mb = [i for i in aliases if i in wes_genes]
            if(len(aliases_mb)>0):
                print("Alias with facets CN value was found for:", k)
                dict_alias[k] = aliases_mb[0]
            else:
                print("No alias with facets CN value was found for:", k)
                missing_value.append(k)
        
        if len(df_missing_genes_no_alias) != 0:
            missing = list(set(list(chain(list(set(df_missing_genes_no_alias)) + list(set(missing_value))))))
        else:
            missing = missing_value
        print("Number of ENiClust genes with no CN values found in FACETS with their names and alias:", len(missing),"/",len(gene_list))
        
        for k, v in dict_alias.items():
            try:
                df_CN[k] = df_CN[v]
            except KeyError:
                pass

        df_CN[missing] = pd.DataFrame(np.nan, index=df_CN.index, columns=missing)
        print("ENiClust genes with at least one sample with no CN value: ",len(df_CN[df_CN.columns.intersection(gene_list)].columns[df_CN[df_CN.columns.intersection(gene_list)].isnull().any()]))
        df_CN = df_CN.fillna(0)
        return(df_CN)

### For 02. segments-related features
    def read_seg_data(self):
        '''This method reads the CCNF file from facets '''
        seg_data = pd.read_csv(self.seg_file, sep = '\t')
        return seg_data

    def get_CN_features(self, hg_build='hg38'):
        '''
        Compute breakpoints per MB by chrom and average: seg_file from facets
        '''
        seg_data = self.read_seg_data()

        ## get hg_build specific chrom sizes
        chrom_sizes = pd.read_csv(constants.hgbuild_map[hg_build], sep = '\t', header = None)

        
        chrom_sizes.columns = ['Chromosome', 'Length']
        chrom_sizes = dict(zip(chrom_sizes['Chromosome'], chrom_sizes['Length']))
        #chrom_sizes = constants.chrom_sizes
        
        if 'ID' in seg_data.columns:
            seg_data['Sample'] = seg_data['ID']
        seg_data = seg_data[['Sample', 'tcn.em', 'chrom', 'loc.start', 'loc.end' ]]
        seg_data['chromosome'] = "chr"+seg_data['chrom'].astype(str)
        seg_data['chrom_length'] = seg_data['chromosome'].map(chrom_sizes)
        seg_data['diff'] = seg_data['tcn.em'].diff()
        seg_data['diff'] =  seg_data['diff'].fillna(0)
        seg_data['CN_change'] = np.where( (seg_data['diff'] >4.0) | (seg_data['diff'] <-4.0), 1, 0) # setting thresholds from seg data
        
        ## obtain  breakpoints per MB per chromosome
        chrom_breakpoints = pd.DataFrame({
            "Num_breakpoints": seg_data.groupby(['Sample', 'chrom', 'chrom_length']).size(),
        }).reset_index()
        chrom_breakpoints['breakpoints_per_mb'] = chrom_breakpoints['Num_breakpoints']/chrom_breakpoints['chrom_length']
        chrom_breakpoints['breakpoints_per_mb'] = chrom_breakpoints['breakpoints_per_mb']*1000000
        
        ##get max breakpoints in chr
        breakpoints_summary = pd.DataFrame({
            "max_breakpoints_any_chr": chrom_breakpoints.groupby('Sample')['breakpoints_per_mb'].max(), # max_breakpoints_any_chr (removed afterwards)
            "avg_breakpoints_chr": chrom_breakpoints.groupby('Sample')['breakpoints_per_mb'].mean(),
        }).reset_index()

        return breakpoints_summary

### For 03. QC features
    def extract_relevant_qc_features(self):
        '''
        Subset the Facets QC file to extract the following info:
        1. Sample
        2. genome_doubled
        3. ploidy
        4. hrd_loh
        5. fraction_cna

        returns: relelavnt_QC_features
        '''
        relevant_QC_features = pd.read_csv(self.facets_qc_file, sep = '\t')
        if 'sample' in relevant_QC_features.columns:
            relevant_QC_features['Sample'] = relevant_QC_features['sample']
        relevant_QC_features = relevant_QC_features[constants.facets_qc_feats]
        return relevant_QC_features

### For 05. receptors status features
    def get_clinical_meta(self):
        '''
        Returns a pd DF for clinical metadata
        Metadata MUST have:
        'Sample', 'ER', "HER2"
        Returns:
            clin_subset: subset of clinical data with Sample, ER and HER2 status
        '''
        clin_info = pd.read_csv(self.clin_meta_file, sep = "\t", encoding = "ISO-8859-1")
        clin_subset = clin_info[['Sample', 'ER', 'HER2']]
        return clin_subset

### Functions to create, scale, transform features
    def _assign_genes_to_ic10_df(self, ic_data):
        '''
        Internal method to correct gene CN levels in case of absence of WGD event 
        and binning of gene CN levels to avoid overfitting
        '''
        ic_data = ic_data.fillna(0)
        
        print("Number of IC1 genes:", len(constants.genes_ic1), "== 63?")
        print("Max total CN value under 500?",  max(ic_data[self.selected_genes['genes'].tolist()].max()) < 500)
        print("Any WGD value missing?", any(pd.isnull(ic_data['genome_doubled'])))
        
        for gene in self.selected_genes['genes']:

            ic_data[gene] = ic_data[gene] * np.where(ic_data['genome_doubled'] == 0, 2, 1)
            
            if gene in constants.genes_ic1: 
                ic_data[gene] = pd.cut(ic_data[gene], constants.bins_IC1, right = False)
            else:
                try:
                    ic_data[gene] = pd.cut(ic_data[gene], constants.bins, right = False)
                except KeyError:
                    continue
            ic_data[gene] = ic_data[gene].astype(str)
            ic_data[gene] = ic_data[gene].str.split(",").str.get(0).str.replace("[", "")
            ic_data[gene] = ic_data[gene].astype(float)
            ic_data[gene] = ic_data[gene].map(constants.dict_bins)

        return ic_data

    def _get_genes_fromCN(self, ic_data, hg_build):
        '''
        Internal method to create the arms features
        '''
        df_ic_10_all = ic_data.copy()
        missing_genes = []
        for chrom_arm in  constants.chrom_cytoband_features: # chromosome arms of interest
            if hg_build == 'hg38':
                genes_in_arm = [i for i in constants.cyto_dict_38.keys() if chrom_arm in constants.cyto_dict_38[i]]
            else:
                genes_in_arm = [i for i in constants.cyto_dict_19.keys() if chrom_arm in constants.cyto_dict_19[i]]
            print("Is it the right version of the ENiClust genes list?", len(list(self.selected_genes['genes'])) == 611)
            genes_in_arm = [i for i in genes_in_arm if i in list(self.selected_genes['genes'])]
            try:
                ic_data[chrom_arm] = ic_data[genes_in_arm].sum(axis=1)
                ic_data[chrom_arm] = ic_data[chrom_arm]/len(genes_in_arm)
                df_ic_10_all[chrom_arm] = ic_data[chrom_arm]

            except KeyError: ##some genes present in chrom_arm level file but don't have copy number data
                genes_with_missing_copycalls = [i for i in genes_in_arm if i not in ic_data.columns.tolist()]
                missing_genes.append(genes_with_missing_copycalls[0])
        return missing_genes,ic_data

    def _transform_ic10_data(self, ic_data,hg_build):
        '''
        Internal method to create, scale, and transform features
        '''
        ic_data = ic_data[~(pd.isnull(ic_data['genome_doubled']))]
        # Scale binary features
        ic_data['genome_doubled'] = 4* ic_data['genome_doubled'].astype(int)
        ic_data['ploidy'] =  ic_data['ploidy']
        ic_data['ER'] = 4* ic_data['ER']
        ic_data['HER2'] = 4* ic_data['HER2']
        
        ic_data['avg_breakpoints_chr'] = 10* ic_data['avg_breakpoints_chr']
        ic_data['fraction_cna'] = 2* ic_data['fraction_cna']
        ic_data['hrd_loh'] = 4* ic_data['hrd_loh']
 
        ic_data = self._assign_genes_to_ic10_df(ic_data) #correction + bins
        missing_genes,ic_data = self._get_genes_fromCN(ic_data,hg_build) #arms
        # ic_data = self._redefine_labels(ic_data, self.keep_clusters, 'ic10')
        ic_data = ic_data.fillna(0)
        return ic_data

### For iC10 reference if required
    def _merge_with_iC10_preds(self, df10):
        '''Internal method
        Returns:
        df10: merged data with iC10 R package predictions
        '''
        dict_ic10_pred = pd.read_csv(self.iC10_results_file, sep = '\t')
        
        dict_ic10_pred = dict_ic10_pred.set_index('Sample', drop = True) 
        dict_ic10_pred = dict_ic10_pred.loc[list(df10['Sample']),:]
        dict_ic10_pred = dict_ic10_pred.reindex(index=df10['Sample'])
        dict_ic10_pred = dict_ic10_pred.reset_index()
        dict_ic10_pred.index = dict_ic10_pred['Sample']
        dict_ic10_pred = dict_ic10_pred[['ic10']]
        
        dict_ic10_pred.columns = ['ic10_exp_CN']
        
        dict_ic10_pred['ic10_exp_CN'] = 'ic'+dict_ic10_pred['ic10_exp_CN'].astype(str)
        dict_ic10_pred = dict(zip(dict_ic10_pred.index, dict_ic10_pred['ic10_exp_CN']))
        
        df10['ic10'] = df10.index.to_series().map(dict_ic10_pred)
        
        df10 = df10[~(pd.isnull(df10['ic10']))]
        df10 = df10.fillna(0)
        return df10

    def _redefine_labels(self, df10, keep_clusters = ['ic1', 'ic2', 'ic6', 'ic9', 'ic5', 'ic3',  'ic8',  'ic10', 'ic4'],
                         ic_field='ic10'):
        df10[ic_field] = np.where(df10[ic_field].isin([ 'ic4ER+', 'ic4ER-']), 'ic4', df10[ic_field])
        df10[ic_field] = np.where(df10[ic_field].isin(keep_clusters), df10[ic_field], "Other")
        labels = df10.pop(ic_field)
        return df10, labels

### For data set split into training(+testing)/validation if required
    def split_training_set(self, ic_data, labels, n_splits = 1, test_size = 0.2):
            '''
            Splits the training set into 80:20
            80: training/test set
            20: validation set
            The training set is used for training the classifier and can be split into training/testing if required
            The validation set is used to validate the classifier and print out the metrics
            '''
            skt = StratifiedShuffleSplit(n_splits = n_splits, test_size = test_size, random_state = 6)

            for train_index, test_index in skt.split(ic_data, labels):
                X_train, X_test = ic_data.iloc[train_index], ic_data.iloc[test_index]
                y_train, y_test = labels.iloc[train_index], labels.iloc[test_index]
            return X_train,X_test,y_train,y_test

### Preprocessing
    def pre_process(self, ic10_results_file=None, split_training_set=False, hg_build='hg38', wgs=False, out_dir="./eniclust_results/preproc"):
        
        gene_CN = self.gene_level_CN()
        gene_CN_subset = self.get_gene_alias(gene_CN, genes_ref_file=constants.training_gene_sets['ALIAS'])
        gene_CN_subset = gene_CN_subset[self.selected_genes['genes']]

        breakpoints_data = self.get_CN_features(hg_build)
        rel_feats = self.extract_relevant_qc_features()
        clin_data = self.get_clinical_meta()

        rel_feats = pd.merge(rel_feats, clin_data, on = 'Sample', how = 'left')
        rel_feats = pd.merge(rel_feats, breakpoints_data, on = 'Sample', how = 'left')
        
        for gene in self.selected_genes['genes'].tolist():
            try:
                gene_CN_subset[gene] = gene_CN[gene]
            except KeyError:
                continue
        gene_CN_subset = gene_CN_subset.drop_duplicates()
        gene_CN_subset = gene_CN_subset[gene_CN_subset.index.isin(rel_feats['Sample'].unique().tolist())]
        
        ic10_data = gene_CN_subset.copy()
        ic10_data['Sample'] = ic10_data.index
        
        ic10_data.index.name = None
        ic10_data = pd.merge(ic10_data, rel_feats, on = 'Sample', how = 'left').set_index(ic10_data['Sample'])
        ic10_data = ic10_data.drop_duplicates()
        
        ic10_data = self._transform_ic10_data(ic10_data, hg_build)

        if ic10_results_file is not None:
            ic10_data = self._merge_with_iC10_preds(ic10_data)
            ic10_data, labels = self._redefine_labels(ic10_data, self.keep_clusters, 'ic10')
         
        if wgs == True:
            ic10_data['WGS'] = 1
            ic10_data['WES'] = 0
        else:
            ic10_data['WGS'] = 0
            ic10_data['WES'] = 1

        ic10_data['Sample'] = ic10_data.index

        if ~(split_training_set):
            
            # Bimodal features
            columns_to_check = ['genome_doubled', 'ER', 'HER2', 'WGS', 'WES']
            print("Are they bimodal? ",all([ic10_data[x].isin([0,1,4]).all() for x in columns_to_check]))
            from collections import Counter
            print("WGD format? ",Counter(ic10_data['genome_doubled']))

        from itertools import chain
        genes_trans = self.selected_genes['genes'].tolist()
        arms_trans = ['11q13', '8q24', '8p11', '17q24']
        other_features_trans = ['genome_doubled', 'ploidy','hrd_loh', 'fraction_cna', 'avg_breakpoints_chr','ER', 'HER2', 'WGS', 'WES']
                    
        # Perform a 80:20 split into training(+testing)/validation data sets
        if split_training_set:
            print('split into training(+testing)/validation')
            
            tcga_samples = [x for x in ic10_data['Sample'] if "TCGA" in x]
            hmf_samples = [x for x in ic10_data['Sample'] if "TCGA" not in x]
            ic10_data.index = ic10_data['Sample']

            cohort_tmp = list()
            cohort_tmp.insert(1, ic10_data[(ic10_data.index.isin(tcga_samples))])
            cohort_tmp.insert(1, ic10_data[(ic10_data.index.isin(hmf_samples))])

            train_ic_10_merged_X = pd.DataFrame()
            test_ic_10_merged_X = pd.DataFrame()
            train_ic_10_merged_y = pd.DataFrame()
            test_ic_10_merged_y = pd.DataFrame()
            
            for df_ic_10 in cohort_tmp:

                print("Sub-cohort size: " + str(df_ic_10.shape[0]))
                
                mysamples = df_ic_10['Sample']
                labels = pd.read_csv(os.path.join(get_data_dir(), 'input/ic10_references.txt'), sep = '\t')
                labels.index = labels['Sample']

                labels = labels[(labels.index.isin(mysamples))]
                labels['ic10'] = np.where(labels['ic10'].isin([ 'ic4ER+', 'ic4ER-' ]), 'ic4', labels['ic10'])
                labels['ic10'] = np.where(labels['ic10'].isin(self.keep_clusters), labels['ic10'], "Other")

                labels = labels.pop('ic10')

                training_data_set,test_data_set,training_set_labels,test_set_labels = self.split_training_set(df_ic_10, labels, n_splits = 1, test_size = 0.2)

                train_ic_10_merged_X = pd.concat([train_ic_10_merged_X, training_data_set ])
                test_ic_10_merged_X = pd.concat([test_ic_10_merged_X, test_data_set ])
                train_ic_10_merged_y = pd.concat([train_ic_10_merged_y, training_set_labels ])
                test_ic_10_merged_y = pd.concat([test_ic_10_merged_y, test_set_labels ])

            training_data_set = train_ic_10_merged_X
            test_data_set = test_ic_10_merged_X
            training_set_labels = train_ic_10_merged_y
            test_set_labels = test_ic_10_merged_y
            
            training_data_set['Sample'] = training_data_set.index
            training_data_set['WGS'] = 0
            training_data_set['WES'] = 0
            training_data_set.loc[[x for x in training_data_set['Sample'] if "TCGA" in x],"WES"] = 1
            training_data_set.loc[[x for x in training_data_set['Sample'] if "TCGA" not in x],"WGS"] = 1
            training_data_set = training_data_set[list(chain(genes_trans + arms_trans + other_features_trans))]
            test_data_set = test_data_set[list(chain(genes_trans + arms_trans + other_features_trans))]
            
            training_data_set = training_data_set.drop(['Sample', 'purity', 'lst', 'ntai', 'max_breakpoints_any_chr'], axis = 1, errors='ignore')

            test_data_set['Sample'] = test_data_set.index
            test_data_set['WGS'] = 0
            test_data_set['WES'] = 0
            test_data_set.loc[[x for x in test_data_set['Sample'] if "TCGA" in x],"WES"] = 1
            test_data_set.loc[[x for x in test_data_set['Sample'] if "TCGA" not in x],"WGS"] = 1
            test_data_set = test_data_set.drop(['Sample', 'purity', 'lst', 'ntai', 'max_breakpoints_any_chr'], axis = 1, errors='ignore')

            print("Train/Test set size:" , training_data_set.shape[0], "samples with" , training_data_set.shape[1], "features")
            print("Validation set size:" , test_data_set.shape[0], "samples with" , test_data_set.shape[1], "features")

            print("Column order maintained:", all(training_data_set.columns == list(chain(self.genes_trans + self.arms_trans + self.other_features_trans))))
            print("Column order maintained:", all(test_data_set.columns == list(chain(self.genes_trans + self.arms_trans + self.other_features_trans))))
            
            training_data_set.to_csv(os.path.join(out_dir,'training_and_testing_set.txt'), sep = '\t', index = True)
            test_data_set.to_csv(os.path.join(out_dir,'validation_set.txt'), sep = '\t', index = True)
            training_set_labels.to_csv(os.path.join(out_dir,'training_and_testing_labels.txt'), sep = '\t', index = True)
            test_set_labels.to_csv(os.path.join(out_dir,'validation_labels.txt'), sep = '\t', index = True)
            return training_data_set, test_data_set, training_set_labels, test_set_labels
        else:
            ic10_data = ic10_data.drop(['Sample', 'purity', 'lst', 'ntai', 'max_breakpoints_any_chr'], axis = 1, errors='ignore')
            ic10_data = ic10_data[list(chain(genes_trans + arms_trans + other_features_trans))]
            print("Column order maintained:", all(ic10_data.columns == list(chain(self.genes_trans + self.arms_trans + self.other_features_trans))))
            ic10_data.to_csv(os.path.join(out_dir,'data_set.txt'), sep = '\t', index = True)
            return ic10_data

def main():
    ec_pp_parser = argparse.ArgumentParser(
        description = 'Module to run the pre-processing steps for ENiClust')

    ec_pp_parser.add_argument('-g', '--gene_copy_number_file',
        help='gene CN level file from Facets [REQUIRED]', 
        action='store', required = True, type=str)
            
    ec_pp_parser.add_argument('-s', '--seg_file',
        help='segments file from Facets [REQUIRED]', 
        action='store', required = True, type=str)

    ec_pp_parser.add_argument('-q', '--facets_qc_file',
        help='Facets QC file from Facets [REQUIRED]', 
        action='store', required = True, type=str)

    ec_pp_parser.add_argument('-f', '--clin_meta_file',
        help='clinical metadata file with receptors status [REQUIRED]', 
        action='store', required = True, type=str)
    
    ec_pp_parser.add_argument('-r', '--ic10_results_file',
        help='iC10 R package predictions file [OPTIONAL]', 
        action='store', required = False, type=str)

    ec_pp_parser.add_argument('-a', '--assembly',
        help='Reference genome [gh19, or hg38]',
        action='store', required = False, type=str, default = 'hg38')

    ec_pp_parser.add_argument('-w', '--wgs',
        help='Set if data is whole genome sequencing',
        action='store_true', default = False)

    ec_pp_parser.add_argument('-y', '--cytoband_map',
        help='exome cytband information [REQUIRED]', 
        action='store', required = False, type=str, default = 'general/exome_gene_locations_arminfo.txt')

    ec_pp_parser.add_argument('-c', '--ic_cluster',
        help='iC10 clusters as classes to predict [OPTIONAL]', 
        action='store', required = False, type=str, default = 'ic1,ic10,ic2,ic4,ic5,ic6,ic8,ic9')

    ec_pp_parser.add_argument('-p', '--split_train',
        help='flag to split the data set into training(+testing)/validation sets [OPTIONAL]', 
        action='store_true')

    ec_pp_parser.add_argument('-o', '--out_dir',
        help='output directory [OPTIONAL]', 
        action='store', required = False, type=str, default = './eniclust_results/preproc')

    args = ec_pp_parser.parse_args()

    gene_copy_number_file = args.gene_copy_number_file
    seg_file = args.seg_file
    facets_qc_file = args.facets_qc_file
    clin_meta_file = args.clin_meta_file
    if args.ic10_results_file:
        ic10_results_file = args.ic10_results_file
    else:
        ic10_results_file = None
    cytoband_map = args.cytoband_map
    keep_clusters = args.ic_cluster.split(',')
    split_training_set = args.split_train

    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        create_dir(out_dir)

    hg_build = args.assembly


    pp = PreProc(gene_copy_number_file, 
                    seg_file, 
                    facets_qc_file, 
                    clin_meta_file,
                    ic10_results_file, 
                    cytoband_map, 
                    keep_clusters)
    
    if split_training_set:
        training_data_set, test_data_set, training_set_labels, test_set_labels = pp.pre_process(ic10_results_file, split_training_set, hg_build, args.wgs, out_dir)
    else:
        processed_data_set = pp.pre_process(ic10_results_file, split_training_set, hg_build, args.wgs, out_dir)

    return None

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("I got interrupted. :-( Bye!\n")
        sys.exit(0)
