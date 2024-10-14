#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on Sep 14, 2022

@author: prisnir, lisemangiante
'''

import os
import pandas as pd
import argparse

import PreProc as pp
import ClassifyIC as ic

from eniclust.utils import get_data_dir, create_dir


### Preprocessing
def format_test_datasset(gene_CN_file, seg_file, facets_qc_file, clin_meta_file, ic10_results_file,
    cytoband_map = os.path.join(get_data_dir(),'general/exome_gene_locations_arminfo.txt'), 
    keep_clusters = ['ic1', 'ic10', 'ic2', 'ic4', 'ic5', 'ic6', 'ic8', 'ic9'], 
    split_training_set = False, hg_build='hg38', wgs=False, out_dir = './eniclust_results'):
    
    tp = pp.PreProc(gene_CN_file, seg_file, facets_qc_file, clin_meta_file, ic10_results_file, cytoband_map, keep_clusters)
    
    if split_training_set:
        training_data_set, test_data_set, training_set_labels, test_set_labels = tp.pre_process(ic10_results_file, split_training_set, hg_build, wgs, out_dir)
        return training_data_set, test_data_set, training_set_labels, test_set_labels
    else:
        processed_data_set = tp.pre_process(ic10_results_file, split_training_set, hg_build, wgs, out_dir)
        return processed_data_set
    
### Running the classifier
def preds(data_set, validation_set_labels = None, 
            keep_clusters = ['ic1', 'ic10', 'ic2', 'ic4', 'ic5', 'ic6', 'ic8', 'ic9'], 
            validate = False,  model_dir = None, out_dir = "./eniclust_results"):

    ec = ic.ClassifyIC(keep_clusters)
    
    if validate:
        preds, report = ec.eniclust_1_0(data_set, validation_set_labels, validate, model_dir, out_dir)
        return preds,report
    else:
        preds = ec.eniclust_1_0(data_set, None, validate, model_dir, out_dir)
        return preds
    
def main():
    # User inputs
    # Overwrite training data defaults with user input
    ec_parser = argparse.ArgumentParser(
        description='ENiClust: Ensemble Integrative Clustering')

    ec_parser.add_argument('-g', '--gene_copy_number_file',
        help='gene CN level file from Facets [REQUIRED]', 
        action='store', required = True, type=str)
            
    ec_parser.add_argument('-s', '--seg_file',
        help='segments file from Facets [REQUIRED]', 
        action='store', required = True, type=str)

    ec_parser.add_argument('-q', '--facets_qc_file',
        help='Facets QC file from Facets [REQUIRED]', 
        action='store', required = True, type=str)

    ec_parser.add_argument('-f', '--clin_meta_file',
        help='clinical metadata file with receptors status [REQUIRED]', 
        action='store', required = True, type=str)
    
    ec_parser.add_argument('-r', '--ic10_results_file',
        help='iC10 R package predictions file [OPTIONAL]', 
        action='store', required = False, type=str)

    ec_parser.add_argument('-y', '--cytoband_map',
        help='exome cytband information [REQUIRED]', 
        action='store', required = False, type=str, default = 'general/exome_gene_locations_arminfo.txt')
    
    ec_parser.add_argument('-c', '--ic_cluster',
        help='iC10 clusters as classes to predict [OPTIONAL]', 
        action='store', required = False, type=str, default = 'ic1,ic10,ic2,ic4,ic5,ic6,ic8,ic9')

    ec_parser.add_argument('-p', '--split_train',
        help='flag to split the data set into training(+testing)/validation first and then, training(+testing) set into training/testing sets [OPTIONAL]', 
        action='store_true')

    ec_parser.add_argument('-v', '--validate',
        help='flag to run the classifier on the validation data set and save a classification report [OPTIONAL]', 
        action='store_true')

    ec_parser.add_argument('-i', '--model_dir', 
        help = 'directory that contains the trained model [OPTIONAL]', 
        action='store', required = False, type=str)

    ec_parser.add_argument('-a', '--assembly',
        help='Reference genome [gh19, or hg38]',
        action='store', required = False, type=str, default = 'hg38')

    ec_parser.add_argument('-w', '--wgs',
        help='Set if data is whole genome sequencing',
        action='store_true', default = False)

    ec_parser.add_argument('-o', '--out_dir',
        help='output directory [OPTIONAL]', 
        action='store', required = False, type=str, default = './eniclust_results')
        
    args = ec_parser.parse_args()
    
    gene_CN_file = args.gene_copy_number_file
    seg_file = args.seg_file
    facets_qc_file = args.facets_qc_file
    clin_meta_file = args.clin_meta_file
    cytoband_map = args.cytoband_map
    ic10_results_file = args.ic10_results_file
    keep_clusters = args.ic_cluster.split(',')
    split_training_set = args.split_train
    validate = args.validate
    hg_build = args.assembly

    if args.model_dir:
        model_dir = os.path.join(get_data_dir(), args.model_dir)
    else:
        model_dir = None

    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        create_dir(out_dir)

    if split_training_set:
        training_data_set, data_set, training_set_labels, validation_set_labels = format_test_datasset(
            gene_CN_file, 
            seg_file, 
            facets_qc_file, 
            clin_meta_file, 
            ic10_results_file, 
            cytoband_map, 
            keep_clusters, 
            split_training_set,
            hg_build,
            args.wgs,
            out_dir)
    else:
        data_set = format_test_datasset(gene_CN_file, seg_file, facets_qc_file, clin_meta_file, ic10_results_file, cytoband_map, keep_clusters, split_training_set, hg_build, args.wgs, out_dir)
    
    if validate:
        ic_preds, report, preds_test, report_test = preds(data_set, validation_set_labels, keep_clusters, validate, model_dir, out_dir)
        pd.DataFrame(preds_test).to_csv(os.path.join(out_dir, 'testing_predictions.txt'), sep = '\t', index = False)
        pd.DataFrame(ic_preds).to_csv(os.path.join(out_dir, 'validation_predictions.txt'), sep = '\t', index = False)
        pd.DataFrame(report_test).to_csv(os.path.join(out_dir, 'classification_report.txt'), sep = '\t', index = True)
        pd.DataFrame(report).to_csv(os.path.join(out_dir, 'classification_report_validation.txt'), sep = '\t', index = True)
    else:
        ic_preds = preds(data_set, None, keep_clusters, validate, model_dir, out_dir)
        pd.DataFrame(ic_preds).to_csv(os.path.join(out_dir, 'data_predictions.txt'), sep = '\t', index = False)
       

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("I got interrupted. :-( Bye!\n")
        sys.exit(0)
