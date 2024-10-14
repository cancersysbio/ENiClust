from collections import defaultdict
import os
import pandas as pd

from eniclust.utils import get_data_dir

RAND = 889

cytoband_map = os.path.join(get_data_dir(), 'general/exome_gene_locations_arminfo.txt')
cyto_info = pd.read_csv(cytoband_map, sep = '\t')

cyto_info['Cytoband_hg19'] = cyto_info['Cytoband_hg19'].str.split(".").str.get(0)
cyto_dict_19 = dict(zip(cyto_info['Gene_symbol'], cyto_info['Cytoband_hg19']))

cyto_info['Cytoband_hg38'] = cyto_info['Cytoband_hg38'].str.split(".").str.get(0)
cyto_dict_38 = dict(zip(cyto_info['Gene_symbol'], cyto_info['Cytoband_hg38']))

selected_genes = pd.read_csv(os.path.join(get_data_dir(), 'general/selected_genes.txt'), sep = '\t')
genes_ic1 = list(selected_genes[selected_genes['ic'] == "ic1"]['genes'])

training_gene_sets = {'EVERYTHING': os.path.join(get_data_dir(), 'general/data_input_genes.txt'),
                      'WES': os.path.join(get_data_dir(), 'general/selected_genes.txt'),
                      'ALIAS': os.path.join(get_data_dir(), 'general/exome_gene_locations_arminfo.txt')
                     }

chrom_cytoband_features = ['11q13', '8q24', '8p11', '17q24']

facets_qc_feats = ['Sample', 'genome_doubled', 'ploidy', 'hrd_loh', 'fraction_cna']

# Bins to transform gene-level CN data
bins = [float(i) for i in [0, 6, 8, 12, 15, 20, 60, 500]]
dict_bins = {k: v for v, k in enumerate(bins)}
dict_bins[4] = 1
bins_IC1 = [float(i) for i in [0, 4, 8, 12, 15, 20, 60, 500]]

hgbuild_map = {
            'hg19':  os.path.join(get_data_dir(), 'general/hg19.chrom.sizes.txt'),
            'hg38':  os.path.join(get_data_dir(), 'general/hg38.chrom.sizes.txt')
            }
