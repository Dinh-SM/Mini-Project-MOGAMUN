# Imports
import glob
import os.path
import moganum_fun
import numpy as np
import pandas as pd
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
		if not os.path.exists(NodesScoresPath):
			nodes_scores = 