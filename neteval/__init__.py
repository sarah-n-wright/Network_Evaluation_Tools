# -*- coding: utf-8 -*-

"""Top-level package for Network Evaluations Tools v2."""

__author__ = """Sarah Wright"""
__email__ = 'snwright@ucsd.edu'
__version__ = '0.2.2'

from .data_import_export_tools import load_public_network_from_ndex, load_private_network_from_ndex, create_networkx_for_export, create_pcnet_networkx_for_export, export_networkx_to_ndex, load_edgelist_to_networkx, write_networkx_to_file, load_node_sets, write_node_sets, 
from .gene_mapper import update_nodes, convert_node_ids, query_mygene
from .query_disgen_associations import create_disgenet_genesets, get_latest_disgenet_disease_list, query_disgenet_disease, query_disgenet, get_disgenet_associations
from .query_gwas_associations import download_file, clean_gwas_catalog_data, create_gwas_gene_sets, process_gwas_genes
from .gsea_functions import GOData, process_go_slim_terms, enrichment_results_from_nodefile, enrich_nodefiles
from .network_evaluation_functions import construct_prop_kernel
from .network_propagation import closed_form_network_propagation, normalize_network
from .network_statistics import NetworkStats, load_network_names
from .node_annotation import load_nodes, get_ensembl_annotation_data, get_ncbi_citation_data, get_uniprot_annotation_data, load_citations, load_ensembl, load_uniprot, load_ensembl, load_hgnc, permutation_test_wrapper, ExpressionData, classify_genes, get_stats_tissue_enriched
from .gene_set_recovery_results import EvaluationResults
from .shuffle_networks import shuffle_network, write_shuffled_network, load_shuffled_network
from .Timer import Timer
from .edge_prediction import EdgePredictionResults, create_network_folds
from .community_annotation import corum_analysis_wrapper, calculate_clustering_coefficient, calculate_homogeneity

import logging

logging.basicConfig(level=logging.INFO)

# Note that import statements could be changed in notebooks for these
# TODO the node annotation functions may require package data. It will be more useful if all annotation data is already included. 

