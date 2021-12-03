# Imports
import mogamun

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

print("START\n")

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

print("DEG, nodes score, layers paths:", deg_path, nodes_scores_path, layers_path, "\n")


evolution_parameters = mogamun.mogamun_init(generations = 10, pop_size = 10, min_size = 10, max_size = 20, crossover_rate = 1, mutation_rate = 1, max_number_of_attemps = 1)

print("Evolution parameters:", evolution_parameters, "\n")

loaded_data = mogamun.mogamun_load_data(evolution_parameters = evolution_parameters, differential_expression_path = deg_path, nodes_scores_path = nodes_scores_path, network_layers_dir = layers_path, layers = "23")

print("Loaded data:", loaded_data, "\n")

random.seed(123)

try:
	res_dir = os.path.abspath("SampleResults/")
except:
	res_dir = None

print("Result dir:", res_dir, "\n")

mogamun.mogamun_run(loaded_data = loaded_data, cores = 2, results_dir = res_dir)