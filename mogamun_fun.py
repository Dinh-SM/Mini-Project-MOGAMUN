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
	merged = pd.DataFrame(list(itertools.zip_longest(V1, V2, Layer)), columns = ["V1", "V2", "Layer"])
	merged = merged.astype({"V1" : str, "V2" : str, "Layer" : str})

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

	search_space_genes = my_mx[0].vs["name"]

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


# Performs a RANDOM depth first search of a given length on a MULTIPLEX network
def make_dfs(
		my_mx,
		root,
		size_of_individual):
	
	# initialization
	discovered = []
	result = []
	current_layer = None
	stack = [root] # initialize the stack of the pending nodes to visit

	# loop while there are pending elements and the desired size hasn't been met
	while (len(root) > 0) and (len(result) < size_of_individual):
		# check if we will change of layer (first loop or 50% chance later)
		if (not current_layer) or (np.random.uniform(0, 1) <= 0.5):
			av_layers = []
			av_layers_ = list(range(1, len(my_mx)+1)) # list all layers
			if len(av_layers_) > 1:
				for av_layer_ in av_layers_:
					if av_layer_ != current_layer:
						av_layers.append(av_layer_)
			else:
				av_layers = av_layers_

			# choose a new layer
			if len(av_layers) > 0:
				current_layer = random.choice(av_layers)
			else:
				current_layer = current_layer

		node = stack[-1] # take a node from the stack to visit it
		stack = stack[:-1] # remove node from stack

		if node not in discovered: # verify if the node hasn't been visited
			discovered.append(node) # add node to the list of visited

			# when finding a disconnected node in a layer, change layer
			keep_going = True
			visited_layers = []
			av_layers = list(range(1, len(my_mx)+1))

			while keep_going:
				my_neighbors_ = my_mx.vs[my_mx.neighbors(node)]["name"]
				my_neighbors = []
				for neighbor in my_neighbors_:
					if neighbor not in discovered:
						my_neighbors.append(neighbor)

				if len(my_neighbors) == 0:
					visited_layers = visited_layers + [current_layer]
					av_layers_ = av_layers + []
					av_layers = []
					for av_layer_ in av_layers_:
						if av_layer_ not in visited_layers:
							av_layers.append(av_layer_)
					if len(av_layers) > 0:
						current_layer = random.choice(av_layers)
					else:
						keep_going = False
				else:
					keep_going = False

			random.shuffle(my_neighbors) # shuffle
			stack = stack + my_neighbors # add shuffled neighbors to stack
			result = result + [node] # add node to the individual

	return result


# Gets the list of IDs of a set of nodes
def get_id_of_nodes(
		list_of_nodes,
		global_network):
	
	nodes = global_network.vs["name"]
	nodes_ids = []
	for node in list_of_nodes:
		nodes_ids.append(nodes.index(node))

	return nodes_ids


# Performs a RANDOM depth first search of a given length on a MULTIPLEX network
def dfs_iterative_mx(
		my_mx,
		root,
		size_of_individual,
		loaded_data):

	max_number_of_attempts = loaded_data["max_number_of_attemps"]
	multiplex = loaded_data["multiplex"]
	keep_looking = True # flag to control the execution of the current function
	attempts = 0

	while keep_looking:
		# make depth first search on MyMx using Root as seed
		result = make_dfs(my_mx, root, size_of_individual)

		# security check of invidual's size
		if len(result) != size_of_individual:
			if attempts > max_number_of_attempts: # if max attemps reached
				my_mx = multiplex # use the big original multiplex
				root = pick_root(my_mx, loaded_data) # pick a random root
			else: # if we still have attempts left
				attempts = attempts + 1 # increment number of attemps
				root = pick_root(my_mx, loaded_data) # pick a root
		else: # if a good root was found and the ind has the desired size
			keep_looking = False # deactivate flag to stop the search

	##### the following conversion has to be done because when a subnetwork is 
	##### created, the nodes get new IDs, according to the size of the new 
	##### subnetwork, but the individuals should always have IDs with respect 
	##### to the global network, therefore the IDs need to be "re-calculated"

	nodes = [] # get local nodes' names
	for res in result:
		nodes.append(my_mx[0].vs["name"][res])

	nodes_ids = get_id_of_nodes(nodes, multiplex[0]) # get IDs

	return nodes_ids


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

		# use DFS (Depth First Search) from the chosen root
		# NOTE. DFS will be used in the initial population because it allows 
		#       to go further from the root

		# makes a random DFS, i.e., the branches are shuffled before
		# visiting them, in order to generate different individuals each time
		dfs = dfs_iterative_mx(multiplex, root, size_of_individual, loaded_data)

		# add individual to the population
		my_population[i] = dfs

	return my_population


# Evaluate an individual
def evaluate_ind(
		individual,
		multiplex,
		loaded_data):

	genes_nodes_scores = loaded_data["genes_with_nodes_scores"]
	density_mx = loaded_data["density_per_layer_multiplex"]

	# verifies that there is at least one node present in the individual
	if len(individual) > 0:

		# gets nodes of the subnetwork
		nodes_of_subnetwork = []
		for i in individual:
			nodes_of_subnetwork.append(multiplex[0].vs["name"][i])

		# gets the sum of the nodes scores of the subnetwork
		genes_nodes_scores_in_subnetwork = genes_nodes_scores[genes_nodes_scores["gene"].isin(nodes_of_subnetwork)]
		genes_nodes_scores_in_subnetwork = list(deg_avail["nodescore"])
		sum_nodes_scores = sum(genes_nodes_scores_in_subnetwork)

		average_nodes_score = sum_nodes_scores / len(individual)

		sum_density_all_layers = 0

		# loop through all the layers in the multiplex
		for layer in range(len(multiplex)):
			# get the inds subnetwork corresponding in the current layer
			subnetwork = multiplex[layer].induced_subgraph(individual)

			#calculate the density of the subnetwork in the current layer
			subnetwork_density = subnetwork.density()

			if not subnetwork_density:
				subnetwork_density = 0

			# add the normalized subnetwork's density with respect to the density of the current layer
			sum_density_all_layers = sum_density_all_layers + (subnetwork_density / density_mx[layer])
	
		res = pd.DataFrame([[average_nodes_score, sum_density_all_layers]], columns = ["average_nodes_score", "density"])
		res = res.astype({"average_nodes_score" : np.float64, "density" : np.float64})
	else:
		res = pd.DataFrame([[0, 0]], columns = ["average_nodes_score", "density"])
		res = res.astype({"average_nodes_score" : np.float64, "density" : np.float64})

	return res


# Evaluates a whole population
def evaluate_population(
		my_population,
		multiplex,
		loaded_data):

	fitness_data = []
	for i in range(len(my_population)):
		fitness_data.append(evaluate_ind(my_population[i], multiplex, loaded_data))

	fitness_data = pd.concat(fitness_data).reset_index(drop = True)

	return fitness_data


# 
def non_dom_sort(
		population,
		loaded_data):

	objective_names = loaded_data["objective_names"]

	#sort individuals by non domination
	ranking = fast_non_dominated_sorting(population[objective_names])

	#TODO


# Defines the function of the body of MOGAMUN
def mogamun_body(
		run_number,
		loaded_data,
		best_inds_path):

	best_inds_file = best_inds_path + "_Run_" + run_number + ".txt"
	my_init_pop = generate_initial_pop(loaded_data["pop_size"], loaded_data["multiplex"], loaded_data)
	fitness_data = evaluate_population(my_init_pop, loaded_data["multiplex"], loaded_data)
	population = pd.DataFrame(data = {}, columns = ["individual"], dtype = object)

	for i in range(len(my_init_pop)):
		df.at[i, "individual"] = my_init_pop[i]

	population = pd.concat([population, fitness_data], axis = 1)

	# obtain ranking and crowding distances
	population = non_dom_sort(population, loaded_data)