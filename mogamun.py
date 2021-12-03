# Imports
import mogamun_fun

import os
import glob
import random
import os.path
import itertools
import numpy as np
import pandas as pd
from igraph import *
import multiprocessing as mp
from scipy.stats import norm
from datetime import datetime
from collections import Counter

# Initialize evolution parameters
def mogamun_init(
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
		"objective_names" : ["average_nodes_score", "density"],
		"measure" : measure,
		"threshold_deg" : threshold_deg,
		"max_number_of_attemps" : max_number_of_attemps
	}

	return evolution_parameters


# Load the data to process
def mogamun_load_data(
		evolution_parameters,
		differential_expression_path,
		nodes_scores_path,
		network_layers_dir,
		layers):

	measure = evolution_parameters["measure"] # "FDR" or "PValue"
	threshold_deg = evolution_parameters["threshold_deg"] # threshold for DEG

	# columns : {gene, logFC, logCPM, PValue, FDR}
	de_results = pd.read_csv(differential_expression_path, dtype = {"gene" : str, "logFC" : np.float64, "logCPM" : np.float64, "PValue": np.float64, "FDR" : np.float64}) # load DE
	de_results = mogamun_fun.remove_duplicates_de_results(de_results) # remove dup entries

	deg = de_results[de_results["FDR"] < threshold_deg] # get list of DEG

	# verify existence of log(fold change) and consider it for the DEG
	if "logFC" in deg.columns:
		deg = deg[abs(deg["logFC"]) > 1]

	# read the file names of the networks for the current experiment, network_layers_dir needs to finish with '/'
	files = glob.glob(network_layers_dir + "/[" + layers + "]_*")

	if len(files) < len(layers):
		print("Error! One or more networks are missing from " + network_layers_dir)
	else:
		# if no nodes scores file exists
		# calculate the nodes scores for all the genes in DE analysis results
		if not os.path.exists(nodes_scores_path):
			nodes_scores = mogamun_fun.get_nodes_scores_of_list_of_genes(de_results, list_of_genes,	measure)

			# data frame of genes and scores. NOTE. Genes not in the list have 0
			genes_with_nodes_scores = pd.DataFrame(list(itertools.zip_longest(de_results["gene"], nodes_scores)), columns = ["gene", "nodescore"]).fillna(0)
			genes_with_nodes_scores = genes_with_nodes_scores.astype({"gene" : str, "nodescore" : np.float64})

			genes_with_nodes_scores.to_csv(nodes_scores_path, index = False)
		else:
			genes_with_nodes_scores = pd.read_csv(nodes_scores_path, dtype = {"gene" : str, "nodescore" : np.float64})

		multiplex = mogamun_fun.generate_multiplex_network(files) # make multiplex network
		merged = mogamun_fun.generate_merged_network(files) # make the merged network
		
		density_per_layer_multiplex = []
		for layer in multiplex:
			density_per_layer_multiplex.append(layer.density())

		loaded_data = evolution_parameters
		loaded_data["network_layers_dir"] = network_layers_dir
		loaded_data["layers"] = layers
		loaded_data["de_results"] = de_results
		loaded_data["deg"] = deg
		loaded_data["genes_with_nodes_scores"] = genes_with_nodes_scores
		loaded_data["multiplex"] = multiplex
		loaded_data["density_per_layer_multiplex"] = density_per_layer_multiplex
		loaded_data["merged"] = merged

		return loaded_data


# Run the algorithm with the specified values for the evolution parameters
def mogamun_run(
		loaded_data,
		cores = 1,
		number_of_runs_to_execute = 1,
		results_dir = '.'):
	
	if loaded_data:
		results_path = results_dir + "/Experiment_" + datetime.today().strftime('%Y-%m-%d') + '/'
		os.mkdir(results_path) # create result folder
		best_inds_path = results_path + "MOGAMUN_Results_" # path for res

		# init with number of cores
		pool = mp.Pool(cores)

		runs = list(range(1, number_of_runs_to_execute+1))

		# run mogamun_body in parallel
		results = [pool.apply(mogamun_fun.mogamun_body, args = (i, loaded_data, best_inds_path)) for i in runs]

		# close pool
		pool.close()

	else:
		print("Missing parameter: loaded_data")