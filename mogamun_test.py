# Imports
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


def main():
	print("START: " + datetime.now().strftime("%H:%M:%S") + "\n")

	try:
		deg_path = os.path.abspath("extdata/DE/Sample_DE.csv")
	except:
		deg_path = None

	try:
		nodes_scores_path = os.path.abspath("extdata/DE/Sample_NodesScore.csv")
	except:
		nodes_scores_path = None

	try:
		layers_path = os.path.abspath("extdata/LayersMultiplex/")
	except:
		layers_path = None

	print("DEG, nodes score, layers paths:", deg_path, nodes_scores_path, layers_path, '\n')


	evolution_parameters = mogamun.mogamun_init(generations = 20, pop_size = 100, min_size = 15, max_size = 50, crossover_rate = 0.8, mutation_rate = 0.1, jaccard_similarity_threshold = 30, tournament_size = 2, measure = "FDR", threshold_deg = 0.05, max_number_of_attemps = 1)

	print("Evolution parameters:", evolution_parameters, '\n')

	loaded_data = mogamun.mogamun_load_data(evolution_parameters = evolution_parameters, differential_expression_path = deg_path, nodes_scores_path = nodes_scores_path, network_layers_dir = layers_path, layers = "23")

	print("Loaded data:", loaded_data, '\n')

	# Plot original network
	plot(loaded_data["merged"], target='./plot/merged_network.png')
	for layer in range(len(loaded_data["multiplex"])):
		plot(loaded_data["multiplex"][layer], target='./plot/multiplex_layer_' + str(layer) + '.png')

	try:
		res_dir = os.path.abspath("SampleResults/")
	except:
		res_dir = None

	print("Result dir:", res_dir, '\n')

	if not os.path.exists(res_dir):
		os.mkdir(res_dir) # create result folder

	mogamun.mogamun_run(loaded_data = loaded_data, cores = 4, results_dir = res_dir, number_of_runs_to_execute = 1)

	# Plot obtained rank 1 subnetworks
	run_number = 1
	result_network_files = glob.glob(os.path.join(res_dir, "Experiment_" + datetime.today().strftime('%Y-%m-%d'), "MOGAMUN_Results__Run_*.txt"))
	for result_network_file in result_network_files:
		result_network = pd.read_csv(result_network_file, names = ["genes", "average_nodes_score", "density", "rank", "crowding_distance"], dtype = {"genes" : object, "average_nodes_score" : np.float64, "density" : np.float64, "rank" : np.intc, "crowding_distance" : np.float64})
		result_network_genes = list(result_network[result_network["rank"] == 1]["genes"])
		for i in range(len(result_network_genes)):
			plot(loaded_data["merged"].induced_subgraph(mogamun_fun.get_id_of_nodes(result_network_genes[i].split(" "), loaded_data["multiplex"][0])), target='./plot/result_network_Run_' + str(run_number) + '_' + str(i) + '.png')
		run_number += 1

if __name__ ==  '__main__':
	import mogamun
	import mogamun_fun
	main()