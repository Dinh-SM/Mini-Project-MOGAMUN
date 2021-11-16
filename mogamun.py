# Imports
import pandas as pd
import moganum_fun

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
		max_number_of_attemps = 3)

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
		layers)

	measure = evolution_parameters["measure"]
	threshold_deg = evolution_parameters["threshold_deg"]

	de_results = pd.read_csv(differential_expression_path)
	de_results = remove_duplicates_de_results(de_results)

	deg = de_results[de_results['FDR'] < threshold_deg]

	# verify existence of log(fold change) and consider it for the DEG
	if "logFC" in deg.columns:
		deg = deg[abs(deg['logFC']) > 1]

	# read the file names of the networks for the current experiment
	