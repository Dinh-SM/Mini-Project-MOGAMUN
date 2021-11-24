# Imports
import glob
import os.path
import moganum_fun
import numpy as np
import pandas as pd
from igraph import *
from scipy.stats import norm

# Initialize evolution parameters
def moganum_init(
		generations = 500,
		pop_size = 100,
		min_size = 15,
		max_size = 50,
		crossover_rate = 0.8,
		mutation_rate = 0.1,
		jaccard_similarity_threshold = 30,
		tournament_size = 2,
		measure = "FDR",
		threshold_deg = 0.05,
		max_number_of_attemps = 3):

	evolution_parameters = {
		"generations" : generations,
		"pop_size" : pop_size,
		"min_size" : min_size,
		"max_size" : max_size,
		"crossover_rate" : crossover_rate,
		"mutation_rate" : mutation_rate,
		"jaccard_similarity_threshold" : jaccard_similarity_threshold,
		"tournament_size" : tournament_size,
		"measure" : measure,
		"threshold_deg" : threshold_deg,
		"max_number_of_attemps" : max_number_of_attemps
	}

	return evolution_parameters


# Load the data to process
def moganum_load_data(
		evolution_parameters,
		differential_expression_path,
		nodes_scores_path,
		network_layers_dir,
		layers):

	measure = evolution_parameters["measure"]
	threshold_deg = evolution_parameters["threshold_deg"]

	# columns : {gene, logFC, logCPM, PValue, FDR}
	de_results = pd.read_csv(differential_expression_path, dtype = {"gene" : str, "logFC" : np.float64, "logCPM" : np.float64, "PValue": np.float64, "FDR" : np.float64})
	de_results = remove_duplicates_de_results(de_results)

	deg = de_results[de_results["FDR"] < threshold_deg]

	# verify existence of log(fold change) and consider it for the DEG
	if "logFC" in deg.columns:
		deg = deg[abs(deg["logFC"]) > 1]

	# read the file names of the networks for the current experiment, network_layers_dir needs to finish with '/'
	files = glob.glob(network_layers_dir + "[" + layers + "]_*")

	if len(files) < len(layers):
		print("Error! One or more networks are missing from " + network_layers_dir)
	else:
		# if no nodes scores file exists
        # calculate the nodes scores for all the genes in DE analysis results
		if not os.path.exists(nodes_scores_path):
			nodes_scores = get_nodes_scores_of_list_of_genes(de_results, list_of_genes,	measure)

			# data frame of genes and scores. NOTE. Genes not in the list have 0
			genes_with_nodes_scores = pd.DataFrame(list(zip(de_results["gene"], nodes_scores)), columns = ["gene", "nodescore"], dtype = {"gene" : str, "nodescore" : np.float64}).fillna(0)

			genes_with_nodes_scores.to_csv(nodes_scores_path, index = False)
		else:
			genes_with_nodes_scores = pd.read_csv(nodes_scores_path, dtype = {"gene" : str, "nodescore" : np.float64})

		multiplex = generate_multiplex_network(files) # make multiplex network