#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on Sep 14, 2022

@author: prisnir, lisemangiante
'''

import pandas as pd
import os
import numpy as np
from collections import Counter
from pyexpat import model
import sys
import pickle
from eniclust import constants
import argparse
from eniclust.utils import get_data_dir

import sklearn.neighbors.base
import sys
sys.modules['sklearn.neighbors._base'] = sklearn.neighbors.base

import sklearn.ensemble.bagging
import sys
sys.modules['sklearn.ensemble._bagging'] = sklearn.ensemble.bagging

import sklearn.ensemble.base
import sys
sys.modules['sklearn.ensemble._base'] = sklearn.ensemble.base

import sklearn.ensemble.forest
sys.modules['sklearn.ensemble._forest'] = sklearn.ensemble.forest

import sklearn.utils.testing
sys.modules['sklearn.utils._testing'] = sklearn.utils.testing

import sklearn.metrics.classification
sys.modules['sklearn.metrics._classification'] = sklearn.metrics.classification

from sklearn import datasets, neighbors, linear_model, ensemble
from sklearn.model_selection import cross_val_score, StratifiedShuffleSplit, StratifiedKFold, RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import classification_report, f1_score, precision_score, recall_score, confusion_matrix, auc, precision_recall_curve 
from sklearn.decomposition import PCA

from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import BorderlineSMOTE
from imblearn.over_sampling import SVMSMOTE
from imblearn.over_sampling import ADASYN

from imblearn.pipeline import Pipeline as imbpipeline
from joblib import Parallel, delayed

import warnings
warnings.filterwarnings('ignore')

class ClassifyIC():
    '''
    Reads training(+testing) data set
    Train the classifier
    Predict results on the testing data set if required and save a classification report
    Predict results on the validation data set and save a classification report or 
    Predict results on new data set
    '''

    def __init__(self, keep_clusters = ['ic1', 'ic10', 'ic2', 'ic4', 'ic5', 'ic6', 'ic8', 'ic9']):
        '''
        Constructor
        '''
        self.feats_excl = ['specimen']
   
        self.keep_clusters = keep_clusters
        self.selected_genes = pd.read_csv(constants.training_gene_sets['WES'], sep = '\t')
        
        self.genes_trans = self.selected_genes['genes'].tolist()
        self.arms_trans = ['11q13', '8q24', '8p11', '17q24']
        self.other_features_trans = ['genome_doubled', 'ploidy','hrd_loh', 'fraction_cna', 'avg_breakpoints_chr', 'ER', 'HER2', 'WGS', 'WES']


### Formatting
    def _prep_incoming_data(self, data, labels = None, keep_clusters = ['ic1', 'ic10', 'ic2', 'ic4', 'ic5', 'ic6', 'ic8', 'ic9']):
        '''
        Internal method to prepare the data set and labels for the classifier
        '''
        if 'Sample' in data.columns:
            data = data.set_index('Sample', drop = True)
        if 'PR' in data.columns:
            data = data.drop('PR', axis = 1)
        from itertools import chain
        print("Column order maintained:", all(data.columns == list(chain(self.genes_trans + self.arms_trans + self.other_features_trans))))
        if labels is not None:
            if isinstance(labels, pd.DataFrame):
                if (len(labels.columns.tolist()) > 1):
                    labels = labels.set_index(labels.columns.tolist()[0], drop = True)
                labels.columns = ['ic10']
                labels['ic10'] = np.where(labels['ic10'].isin(keep_clusters),labels['ic10'],'Other')
                labels = labels['ic10']
                labels = pd.Series(labels)
            return data,labels
        else:
            return data

### Loading the Gradient Booster classifier
    def _init_GBT(self, model_dir=None, out_dir = "./eniclust_results/classifyIC", min_samples_leaf= 1, max_depth = 2, learning_rate = 0.2, random_state=constants.RAND):
        
        ''' 
        Internal method to train or load the GBT classifier
        ''' 
        files = os.listdir(model_dir)
        files = [x for x in files if "gbt__" in x][0]
        path = model_dir + "/" + files
        print("GBT", path)
        with open(path, 'rb') as gbt_model:
             gbt_model = pickle.load(gbt_model)
        
        return gbt_model

### Loading the Random Forest classifier
    def _init_RF(self, model_dir = None, out_dir = "./eniclust_results/classifyIC", 
                 n_estimators= 50, min_samples_leaf= 1, 
                 max_depth= 70, random_state=constants.RAND):
        ''' 
        Internal method to train or load the RF classifier
        '''
        files = os.listdir(model_dir)
        files = [x for x in files if "rf__" in x][0]
        path = model_dir + "/" + files
        print("RF:", path)
        with open(path, 'rb') as rf_model:
             rf_model = pickle.load(rf_model)
        
        return rf_model

### Loading the Elastic Network classifier
    def _init_ENET(self, model_dir = None, out_dir = "./eniclust_results/classifyIC", l1_ratio= 0.9, C= 0.2, random_state=constants.RAND):
        
        ''' 
        Internal method to train or load the ENET classifier
        '''
        files = os.listdir(model_dir)
        files = [x for x in files if "enet__" in x][0]
        path = model_dir + "/" + files
        print("ENET:", path)
        with open(path, 'rb') as enet_model:
             enet_model = pickle.load(enet_model)
        
        return enet_model

### Loading the Voting consensus classifier
    def _init_voting_process(self, model_dir = None, out_dir = "data/test/output_classifyIC"):
        '''
        Internal method to train or load the Voting consensus classifier
        '''
        files = os.listdir(model_dir)
        files = [x for x in files if "voting__" in x][0]
        path = model_dir + "/" + files
        print("Voting:", path)
        with open(path, 'rb') as model:
             model = pickle.load(model)
        
        return model

### Running the classifier
    def eniclust_1_0(self, data_set, validation_set_labels = None, validate = False, model_dir = None, out_dir = "./eniclust_results/classifyIC"):

        if validate:
            validation_set, validation_set_labels = self._prep_incoming_data(data_set, validation_set_labels)
        else:
            validation_set = self._prep_incoming_data(data_set, None)

        preds = pd.DataFrame(columns=["Sample"])
        preds['Sample'] = validation_set.index.values

        if validate:
            preds['ic10'] = validation_set_labels.tolist()
            report = pd.DataFrame(columns=['precision', 'recall', 'f1-score', 'support', 'data', 'model'])

        print ('GBT classifier')          
        gbt_class=self._init_GBT(model_dir, out_dir)

        p = gbt_class.predict(validation_set)
        preds['GBT'] = p
        proba = pd.DataFrame(gbt_class.predict_proba(validation_set))
        proba.columns = ['GBT_proba_Other'] + sorted(["GBT_proba_" + s for s in self.keep_clusters])
        preds = pd.concat([preds, proba], axis = 1)

        if validate:
            report_tmp = classification_report(preds['ic10'], preds['GBT'], output_dict=True)
            report_tmp = pd.DataFrame(report_tmp).transpose()
            report_tmp['data'] = "validation"
            report_tmp['model'] = 'GBT'
            report = pd.concat([report, report_tmp])

        print ('RF classifier')
        rf_class=self._init_RF(model_dir, out_dir)

        p = rf_class.predict(validation_set)
        preds['RF']= p
        proba = pd.DataFrame(rf_class.predict_proba(validation_set))
        proba.columns = ['RF_proba_Other'] + sorted(["RF_proba_" + s for s in self.keep_clusters])
        preds = pd.concat([preds, proba], axis = 1)

        if validate:
            report_tmp = classification_report(preds['ic10'], preds['RF'], output_dict=True)
            report_tmp = pd.DataFrame(report_tmp).transpose()
            report_tmp['data'] = "validation"
            report_tmp['model'] = 'RF'
            report = pd.concat([report, report_tmp])

        print ('ENET classifier')
        enet_class = self._init_ENET(model_dir, out_dir)

        p = enet_class.predict(validation_set)
        preds['ENET']= p
        proba = pd.DataFrame(enet_class.predict_proba(validation_set))
        proba.columns = ['ENET_proba_Other'] + sorted(["ENET_proba_" + s for s in self.keep_clusters])
        preds = pd.concat([preds, proba], axis = 1)

        if validate:
            report_tmp = classification_report(preds['ic10'], preds['ENET'], output_dict=True)
            report_tmp = pd.DataFrame(report_tmp).transpose()
            report_tmp['data'] = "validation"
            report_tmp['model'] = 'ENET'
            report = pd.concat([report, report_tmp])

        print ('Voting classifier')
        eni_model = self._init_voting_process(model_dir, out_dir)

        p = eni_model.predict(validation_set)
        preds['VotingClassifier']= p
        proba = pd.DataFrame(eni_model.predict_proba(validation_set))
        proba.columns = ['Voting_proba_Other'] + sorted(["Voting_proba_" + s for s in self.keep_clusters])
        preds = pd.concat([preds, proba], axis = 1)

        if validate:
            report_tmp = classification_report(preds['ic10'], preds['VotingClassifier'], output_dict=True)
            report_tmp = pd.DataFrame(report_tmp).transpose()
            report_tmp['data'] = "validation"
            report_tmp['model'] = 'VotingClassifier'
            report = pd.concat([report, report_tmp])

            print ('predictions and classification report READY')
            return preds,report
        else:
            return preds

def main():
    # Set up some default parameters
    data_dir= os.path.join(get_data_dir(), 'test/output_preproc')
    training_data_file = os.path.join(data_dir, 'training_set.txt')
    training_data_labels = os.path.join(data_dir, 'training_labels.txt')

    ec_ic_parser = argparse.ArgumentParser(
        description = 'Module to train and run the ENiClust classifier on either the validation or user data set')
    
    ec_ic_parser.add_argument('-t', '--training_set', 
        help = 'training data set file [OPTIONAL]', 
        action='store', required = False, type=str, default = training_data_file)
    
    ec_ic_parser.add_argument('-l', '--training_set_labels', 
        help = 'training data labels file [OPTIONAL]', 
        action='store', required = False, type=str, default = training_data_labels)
    
    ec_ic_parser.add_argument('-d', '--data_set', 
        help = 'data set file from which the classifier will predict the clusters, it can be either validation set or any other new data set [REQUIRED]', 
        action='store', required = True, type=str)
    
    ec_ic_parser.add_argument('-c', '--ic_cluster',
        help='iC10 clusters as classes to predict [OPTIONAL]', 
        action='store', required = False, type=str, default = 'ic1,ic10,ic2,ic4,ic5,ic6,ic8,ic9')
    
    ec_ic_parser.add_argument('-v', '--validate',
        help='flag to run the classifier on the validation data set and save a classification report [OPTIONAL]', 
        action='store_true')
    
    ec_ic_parser.add_argument('-i', '--model_dir', 
        help = 'directory that contains the trained model [OPTIONAL]', 
        action='store', required = False, type=str)
    
    ec_ic_parser.add_argument('-o', '--out_dir', 
        help = 'output directory [REQUIRED]', 
        action='store', required = False, type=str, default = './eniclust_results/classifyIC')
    
    args = ec_ic_parser.parse_args()
    
    data_set_file = args.data_set
    data_set = pd.read_csv(data_set_file, sep = '\t')
    
    keep_clusters = args.ic_cluster.split(',')
    validate = args.validate
    if validate:
        validation_set_labels = os.path.join(data_dir, 'validation_labels.txt') 
        validation_set_labels = pd.read_csv(validation_set_labels, sep = '\t')

    if args.model_dir:
        #model_dir = os.path.join("./", args.model_dir)
        model_dir = os.path.join(get_data_dir(), args.model_dir)
        
    else:
        model_dir = None

    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    ec = ClassifyIC(keep_clusters)

    if validate:
        preds, report = ec.eniclust_1_0(data_set, validation_set_labels, validate, model_dir, out_dir)
        pd.DataFrame(preds).to_csv(os.path.join(out_dir, 'validation_predictions.txt'), sep = '\t', index = False)
        pd.DataFrame(report).to_csv(os.path.join(out_dir, 'classification_report_validation.txt'), sep = '\t', index = True)
    else:
        preds = ec.eniclust_1_0(data_set, None, validate, model_dir, out_dir)
        pd.DataFrame(preds).to_csv(os.path.join(out_dir, 'data_predictions.txt'), sep = '\t', index = False)
    return None

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("I got interrupted. :-( Bye!\n")
        sys.exit(0)
