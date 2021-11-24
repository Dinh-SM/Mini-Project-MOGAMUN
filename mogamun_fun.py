# Filters the DE results, to remove duplicates
def remove_duplicates_de_results(
		de_results):

	# remove NAs
	de_results = de_results.dropna(subset = ["gene"])

	# remove duplicates by keeping the one with the minimal FDR
	de_results_filtered = de_results.sort_values("FDR", ascending = True).drop_duplicates("gene").sort_index()

	return de_results_filtered


# Calculates the node score of a list of genes
def get_nodes_scores_of_list_of_genes(
		de_results,
		list_of_genes,
		measure):
	
	list_of_genes = list(set(list_of_genes))

	# calculate nodes scores for the genes with formula: z = inverse CDF(1-p) 
	# if gene has an expression value, and z = -Inf otherwise
	nodes_scores = []
	for gene in list_of_genes:
		if gene in de_results["gene"].tolist():
			de_results_gene = de_results[de_results["gene"] == gene]
			node_score = norm.ppf(1 - de_results_gene[measure])
		else:
			node_score = np.NINF
		nodes_scores.append(node_score)

	# normalizes nodes scores to be in the range [0-1], ignoring +Inf and -Inf
	for node_score in nodes_scores:
		min_nodes_scores = np.nanmin(nodes_scores[nodes_scores != np.NINF])
		max_nodes_scores = np.nanmax(nodes_scores[nodes_scores != np.INF])
		node_score = (node_score - min_nodes_scores) / (max_nodes_scores - min_nodes_scores)

	# replace INF with 1 and replace NINF with 0
	for node_score in nodes_scores:
		if node_score == np.NINF:
			node_score = 0
		if node_score == np.INF:
			node_score = 1

	return nodes_scores


# Generates the multiplex network
def generate_multiplex_network(
		files):

	# declare empty list to store the multiplex
	multiplex = []

	all_nodes_to_take = get_list_of_all_nodes_present_in_layers(files)

	# loop through all the layers to get the corresponding subnetwork
	for layer_file in files:
		# load layer
		layer = pd.read_csv(layer_file, sep = '\t')

		# create network with the current layer
		# NOTE. All the nodes will have the same ID in every layer
		current_network = Graph(directed = False)
		for node in all_nodes_to_take:
			current_network.add_vertex(str(node))
		for index, row in layer.iterrows():
			current_network.add_edge(str(row[0]), str(row[1])) 
		current_network.vs["label"] = current_network.vs["name"]

		# add subnetwork as a layer into the multiplex network
		multiplex.append(current_network)

	return multiplex


# Gets the sorted list of nodes
def get_list_of_all_nodes_present_in_layers(
		files):

	# declare empty variable to store the nodes
	all_nodes = []

	# loop through all the layers to get the corresponding subnetwork
	for layer_file in files:
		# load layer
		layer = pd.read_csv(layer_file, sep = '\t')

		all_nodes = all_nodes + list(set(pd.concat([layer.iloc[:,0], layer.iloc[:,1]], join = "outer")))

	# remove duplicates and sort alphabetically
	all_nodes = sorted(list(set(all_nodes)))

	# cast to string
	all_nodes = [str(i) for i in all_nodes]

	return all_nodes


# Generates the merged network
def generate_merged_network(
		files):

	# declare empty lists for the edges of the merged
	V1 = []
	V2 = []
	Layer = []

	# loop through all the layers to get all the interactions
	for layer_file in files:
		# load layer
		layer = pd.read_csv(layer_file, sep = '\t')

		V1 = V1 + list(layer.iloc[:,0])
		V2 = V2 + list(layer.iloc[:,1])
		Layer = Layer + [layer_file]*len(list(layer.iloc[:,0]))

	# get the data frame containing the lists
	merged = pd.DataFrame(list(zip(V1, V2, Layer)), columns = ["V1", "V2", "Layer"], dtype = {"V1" : str, "V2" : str, "Layer" : str})

	# create merged network, keeping the same order of the nodes as in the mx
	merged_network = Graph(directed = False)
	for index, row in merged.iterrows():
		merged_network.add_vertex(str(row[0]))
		merged_network.add_vertex(str(row[1]))
		merged_network.add_edge(str(row[0]), str(row[1]))
	merged_network.vs["label"] = merged_network.vs["name"]

	return merged_network


#
def pick_root(
		my_mx,
		loaded_data):
	
	deg = loaded_data["deg"]

	search_space_genes = my_mx.vs["name"]

	# get the list of DE genes in the search space
	deg_avail = deg[deg["gene"].isin(search_space_genes)]
	deg_avail = list(deg_avail["gene"])

	# verify if there is at least one DE gene
	if len(deg_avail) > 0:
		# randomly pick a DEG to be the root of the tree (individual)
		root = deg_avail[random.randint(0, len(deg_avail)-1)]
	else:
		# randomly pick any gene from the list of available genes
		root = search_space_genes[random.randint(0, len(search_space_genes)-1)]

	return root


# Generate the initial population
def generate_initial_pop(
		pop_size,
		multiplex,
		loaded_data):

	min_size = loaded_data["min_size"]
	max_size = loaded_data["max_size"]

	my_population = [0]*pop_size
	for i in range(1, pop_size+1):
		size_of_individual = random.randint(min_size, max_size)

		root = pick_root(multiplex, loaded_data)


#
def evaluate_population():
	pass


# Defines the function of the body of MOGAMUN
def mogamun_body(
		run_number,
		loaded_data,
		best_inds_path):

	best_inds_file = best_inds_path + "_Run_" + run_number + ".txt"
	my_init_pop = generate_initial_pop(loaded_data["pop_size"], loaded_data["multiplex"], loaded_data)
	fitness_data = evaluate_population(my_init_pop, loaded_data["multiplex"], loaded_data)