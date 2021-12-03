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
		max_nodes_scores = np.nanmax(nodes_scores[nodes_scores != np.inf])
		node_score = (node_score - min_nodes_scores) / (max_nodes_scores - min_nodes_scores)

	# replace INF with 1 and replace NINF with 0
	for node_score in nodes_scores:
		if node_score == np.NINF:
			node_score = 0
		if node_score == np.inf:
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
	v1 = []
	v2 = []
	layer = []

	# loop through all the layers to get all the interactions
	for layer_file in files:
		# load layer
		layer_ = pd.read_csv(layer_file, sep = '\t')

		v1 = v1 + list(layer_.iloc[:,0])
		v2 = v2 + list(layer_.iloc[:,1])
		layer = layer + [layer_file]*len(list(layer_.iloc[:,0]))

	# get the data frame containing the lists
	merged = pd.DataFrame(list(itertools.zip_longest(v1, v2, layer)), columns = ["v1", "v2", "Layer"])
	merged = merged.astype({"v1" : str, "v2" : str, "layer" : str})

	# create merged network, keeping the same order of the nodes as in the mx
	merged_network = Graph(directed = False)
	for index, row in merged.iterrows():
		merged_network.add_vertex(str(row[0]))
		merged_network.add_vertex(str(row[1]))
		merged_network.add_edge(str(row[0]), str(row[1]))
	merged_network.vs["label"] = merged_network.vs["name"]

	return merged_network


# Picks a random node to be the root of a search. DEG are preferred over non-DEG
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
			av_layers_ = list(range(len(my_mx))) # list all layers
			if len(av_layers_) > 1:
				for av_layer_ in av_layers_:
					if av_layer_ != current_layer:
						av_layers.append(av_layer_)
			else:
				av_layers = av_layers_

			# choose a new layer
			if len(av_layers) > 0:
				current_layer = random.choice(av_layers)

		node = stack[-1] # take a node from the stack to visit it
		stack = stack[:-1] # remove node from stack

		if node not in discovered: # verify if the node hasn't been visited
			discovered.append(node) # add node to the list of visited

			# when finding a disconnected node in a layer, change layer
			keep_going = True
			visited_layers = []
			av_layers = list(range(len(my_mx)))

			while keep_going:
				my_neighbors_ = my_mx[current_layer].vs[my_mx[current_layer].neighbors(node)]["name"]
				my_neighbors = []
				for neighbor in my_neighbors_:
					if neighbor not in discovered:
						my_neighbors.append(neighbor)

				if len(my_neighbors) == 0:
					visited_layers.append(current_layer)
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
			result.append(node) # add node to the individual

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

	my_population = [[]]*pop_size
	for i in range(pop_size):
		size_of_individual = random.randint(min_size, max_size)

		root = pick_root(multiplex, loaded_data)

		# use DFS (Depth First Search) from the chosen root
		# NOTE. DFS will be used in the initial population because it allows 
		# to go further from the root

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


# Performs the fast non dominated sorting 
def fast_non_dominated_sorting(
		input_data):
	
	pop_size = input_data.shape[0]
	idx_dominators = [[]]*pop_size
	idx_dominatees = [[]]*pop_size

	for i in range(pop_size - 1):
		for j in range(i, pop_size):
			if i != j:
				xi = list(input_data.iloc[i])
				xj = list(input_data.iloc[j])

				if_all_1 = True
				if_any_1 = False
				for k, l in zip(xi, xj):
					if k < l:
						if_all_1 = False
					if k > l:
						if_any_1 = True

				if_all_2 = True
				if_any_2 = False
				for k, l in zip(xi, xj):
					if l < k:
						if_all_2 = False
					if l > k:
						if_any_2 = True

				if if_all_1 and if_any_1:
					idx_dominators[j].append(i)
					idx_dominatees[i].append(j)
				elif if_all_2 and if_any_2:
					idx_dominators[i].append(j)
					idx_dominatees[j].append(i)

	no_dominators = [len(i) for i in idx_dominators]
	rnk_list = [[i for i, no in enumerate(no_dominators) if no == 0]]
	sol_assigned = [len([i for i, no in enumerate(no_dominators) if no == 0])]
	while sum(sol_assigned) < pop_size:
		q = []
		no_sol_in_curr_frnt = sol_assigned[-1]
		for i in range(no_sol_in_curr_frnt):
			sol_idx = rnk_list[-1][i]
			his_dominatees = idx_dominatees[sol_idx]
			for j in his_dominatees:
				no_dominators[j] = no_dominators[j] - 1
				if no_dominators[j] == 0:
					q.append(j)

		rnk_list.append(sorted(q))
		sol_assigned.append(len(q))

	return rnk_list


# Estimates the density of solutions surrounding a particular solution within each front. It calculates the crowding distances of solutions according to their objectives and those within the same front.
def crowding_dist_4_frnt(
		population_,
		ranking,
		range_):

	pop_size = len(population_)
	obj_dim = len(range_)
	var_no = len(population_[0]) - 1 - obj_dim
	cd = []
	for i in range(pop_size):
		cd.append([np.inf]*obj_dim)
	for i in range(len(ranking)):
		select_row = []
		for row in population_:
			if row[-1] == i+1:
				select_row.append(row)
		length = len(ranking[i])
		if length > 2:
			for j in range(obj_dim):
				select_ = []
				for row in select_row:
					select_.append(row[var_no+j])
				select_ = sorted(range(len(select_)), key = lambda k: select_[k])
				original_idx = []
				for idx in select_:
					original_idx.append(ranking[i][idx])
				for k in range(2, length):
					cd[original_idx[k-1]][j] = abs(population_[original_idx[k]][var_no+j] - population_[original_idx[k-2]][var_no+j]) / range_[j]

	return cd


# Performs the fast non dominated sorting and the calculus of the crowding distance
def non_dom_sort(
		population,
		loaded_data):

	objective_names = loaded_data["objective_names"]

	#sort individuals by non domination
	ranking = fast_non_dominated_sorting(population[objective_names])

	# transform the output of the sorting into a dataframe of 2 columns:
	# 1.- Individual ID. 2.- Rank
	individuals = []
	ranks = []
	for i in range(len(ranking)):
		for ind in ranking[i]:
			individuals.append(ind)
			ranks.append(i+1)

	my_result = pd.DataFrame(list(itertools.zip_longest(individuals, ranks)), columns = ["individual", "rank"])

	# order the dataframe by individual ID
	my_result = my_result.sort_values(by = 'individual', ignore_index = True)

	# add the rank to the data frame
	population = population.join(my_result["rank"])
	
	# calculate (MAX - MIN) of every objective function
	range_ = []
	for objective in objective_names:
		objective_col = list(population[objective])
		range_.append(max(objective_col) - min(objective_col))

	# create a dataframe removing the ind codes and the crowding distances
	pop_col = [col for col in list(population.columns) if col not in ["individual", "crowding_distance"]]
	pop_ = population[pop_col]
	population_ = []
	for i in range(pop_.shape[0]):
		population_.append(list(pop_.iloc[i]))

	crowding_distances = crowding_dist_4_frnt(population_, ranking, range_)
	crowding_distance = []
	for cd in crowding_distances:
		crowding_distance.append(sum(cd))

	population["crowding_distance"] = crowding_distance

	return population


# Performs a tournament selection
def tournament_selection(
		tournament_size,
		tournament_pop):

	# randomly choose as many individuals as the tournament size indicates
	ids = []
	while len(ids) != tournament_size:
		id_ = random.choice(list(range(tournament_pop.shape[0])))
		if id_ not in ids:
			ids.append(id_)

	# verify all both individuals are in the same Pareto front
	ranks = list(tournament_pop["rank"])
	ranks_ids = []
	for id_ in ids:
		ranks_ids.append(ranks[id_])

	if len(list(set(ranks_ids))) == 1:
		# verify if these individuals have information about crowding distance
		# (they won't if generation == 1)
		if "crowding_distance" in list(tournament_pop.columns):
			# get id of the individual with the highest crowding distance
			crowding_distances = list(tournament_pop["crowding_distance"])
			crowding_distance_ids = []
			for id_ in ids:
				crowding_distance_ids.append(crowding_distances[id_])
			winner = crowding_distances.index(max(crowding_distance_ids))
		else:
			# if there is no crowding distance and the individuals have the
			# same rank, they are all considered as winners
			winner = list(range(tournament_size))
	else:
		# if the individuals are in different Pareto fronts, the winner
		# of the tournament is the one with lowest value (best ranking)
		winner = ranks.index(min(ranks_ids))

	tournament_winner = tournament_pop.iloc[[winner]]

	return tournament_winner


# Gets the neighbors of a list of nodes, considering all the layers
def get_neighbors(
		node_list,
		multiplex):
	
	neighbors = []

	# loop through all the layer of the multiplex
	for layer in range(len(multiplex)):
		# add to the list the neighbors of the node list in the current layer
		for node in node_list:
			neighbors = neighbors + layer.vs[layer.neighbors(node)]["name"]

	# sort result and delete duplicates
	neighbors = sorted(list(set(neighbors)))

	return neighbors


# Chooses a compatible individual with parent 1
def get_parent2(
		parent1,
		population,
		loaded_data):

	multiplex = loaded_data["multiplex"]
	pop_size = loaded_data["pop_size"]
	tournament_size = loaded_data["tournament_size"]

	### filter population to leave only inds near parent 1
	# get nodes ids with respect to the big network
	nodes_ids_of_parent1 = []
	for ind in list(parent1["individual"]):
		nodes_ids_of_parent1 = nodes_ids_of_parent1 + ind
	nodes_ids_of_parent1 = sorted(nodes_ids_of_parent1)

	# get list of nodes IDs in parent 1 and their neighbors
	neighbors_nodes_p1 = get_neighbors(nodes_ids_of_parent1, multiplex)
	
	# get inds that contain at least one node from the previous list
	neighbors_inds_p1 = []
	individuals = []
	for i in range(pop_size):
		if len(list(set(population.at[i, "individual"]) & set(neighbors_nodes_p1))) > 0:
			neighbors_inds_p1.append(i)

	# filter population and leave individuals near parent 1
	potential_inds_parent2 = population.iloc[neighbors_inds_p1]

	# verify if parent 1 is in the list and if so, remove it
	for index, row in potential_inds_parent2.iterrows():
		for index_, row_ in parent1.iterrows():
			if row == row_:
				potential_inds_parent2.drop([index])
				break

	# if there are more than 2 inds
	if potential_inds_parent2.shape[0] >= 2:
		# tournament to choose parent 2
		parent2 = tournament_selection(tournament_size, potential_inds_parent2)
	elif potential_inds_parent2.shape[0] == 1:
		# if there is only 1 ind, keep it
		parent2 = potential_inds_parent2
	else:
		# unsuccessful search for compatible parents
		parent2 = None

	return parent2


# filters the multiplex network, to keep only the specified list of nodes
def filter_multiplex(
		multiplex,
		nodes_to_keep):

	# declare empty list to store the multiplex
	filtered_multiplex = []

	# loop through all the layers to get the corresponding subnetwork
	for layer in range(len(multiplex)):
		# create network with the current layer
		current_network = multiplex[layer].induced_subgraph(nodes_to_keep)

		# add subnetwork as a layer into the multiplex network
		filtered_multiplex.append(current_network)

	return filtered_multiplex


# perform a RANDOM breadth first search of a given length on a MULTIPLEX network. This means that the branch to be visted is randomnly picked from all the available ones
def make_bfs(
		my_mx,
		root,
		size_of_individual):
	
	# initialization
	discovered = []
	result = []
	current_layer = None
	# initialize the queue of the pending nodes to visit
	queue = [root]

	# loop while there are pending elements and the desired size hasn't been met
	while len(queue) > 0 and len(result) < size_of_individual:
		# check if we will change of layer
		if not current_layer or np.random.uniform(0, 1) <= 0.5:
			av_layers = []
			av_layers_ = list(range(len(my_mx))) # list all layers
			if len(av_layers_) > 1:
				for av_layer_ in av_layers_:
					if av_layer_ != current_layer:
						av_layers.append(av_layer_)
			else:
				av_layers = av_layers_

			# choose a new layer
			if len(av_layers) > 0:
				current_layer = random.choice(av_layers)

		# take a node from the queue to visit it
		node = queue[0]
		queue = queue[1:]

		# verify if the node hasn't been visited
		if node not in discovered:
			# add node to the list of visited
			discovered.append(node)

			# when findind a disconnected node in a layer, change layer
			keep_going = True
			visited_layers = []
			av_layers = list(range(len(my_mx)))

			while keep_going:
				my_neighbors = my_mx[current_layer].vs[my_mx[current_layer].neighbors(node)]["name"]
				my_neighbors = []
				for neighbor in my_neighbors_:
					if neighbor not in discovered:
						my_neighbors.append(neighbor)

				if len(my_neighbors) == 0:
					visited_layers.append(current_layer)
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
			queue = queue + my_neighbors # add shuffled neighbors to queue
			result.append(node) # add node to the individual

	return my_result


# Performs a RANDOM breadth first search of a given length on a MULTIPLEX network. This means that the nodes to be visted are always randomnly picked from all the available ones
def bfs_iterative_mx(
		my_mx,
		root,
		size_of_individual,
		loaded_data):
	
	max_number_of_attempts = loaded_data["max_number_of_attemps"]
	multiplex = loaded_data["multiplex"]
	# flag to control the execution of the current function
	keep_looking = True
	attempts = 0

	while keep_looking:
		result = make_bfs(my_mx, root, size_of_individual)

		# security check of individual's size
		if len(result) != size_of_individual:
			# if max attemps reached
			if attempts > max_number_of_attempts:
				# use the big original multiplex
				my_mx = multiplex
				# pick a random root
				root = pick_root(my_mx, loaded_data)
			# if we still have attemps left
			else:
				# increment number of attemps
				attempts = attempts + 1
				# pick a root
				root = pick_root(my_mx, loaded_data)
		# if a good root was found and the ind has the desired size
		else:
			# deactivate flag to stop the search
			keep_looking = False

	##### the following conversion has to be done because when a subnetwork is 
	##### created, the nodes get new IDs, according to the size of the new 
	##### subnetwork, but the individuals should always have IDs with respect 
	##### to the global network, therefore the IDs need to be "re-calculated"

	# get the names of the nodes in the local network with the local IDs
	nodes = []
	for res in result:
		nodes.append(my_mx[0].vs["name"][res])

	# get the global IDs of the corresponding nodes
	nodes_ids = get_id_of_nodes(nodes, multiplex[0])

	return nodes_ids


# Performs crossover
def crossover(
		parent1,
		parent2,
		loaded_data):
	
	max_number_of_attempts = loaded_data["max_number_of_attemps"]
	min_size = loaded_data["min_size"]
	max_size = loaded_data["mac_size"]
	crossover_rate = loaded_data["crossover_rate"]
	multiplex = loaded_data["multiplex"]
	# intialize empty list
	children = []
	# generate a random number between 0 and 1
	p = np.random.uniform(0, 1)

	# check if crossover is to be performed
	if p <= crossover_rate:
		# join both parents
		nodes_p = []
		for ind1 in list(parent1["individual"]):
			nodes_p = nodes_p + ind1
		for ind2 in list(parent2["individual"]):
			for ind2_ in ind2:
				if ind2_ not in nodes_p:
					nodes_p.append(ind2_)

		# first, generate the multiplex of the joint parents
		mx_parents = filter_multiplex(multiplex, nodes_p)
		
		# if length is minimum, copy parents
		if len(nodes_p) == min_size:
			children = [[], []]
			for ind1 in list(parent1["individual"]):
				children[0] = children[0] + ind1
			for ind2 in list(parent2["individual"]):
				children[1] = children[1] + ind2
		else:
			# loop to generate two children
			for k in range(2):
				# pick a size for the child
				if len(nodes_p) >= max_size:
					size = random.randint(min_size, max_size)
				else:
					size = random.randint(min_size, len(nodes_p))

				# pick root
				root = pick_root(mx_parents, loaded_data)
				# pick search method
				search_method = random.choice(["DFS", "BFS"])

				# depth first search
				if search_method == "DFS":
					dfs = dfs_iterative_mx(mx_parents, root, size, loaded_data)
				# breadth first search
				elif search_method == "BFS":
					bfs = bfs_iterative_mx(mx_parents, root, size, loaded_data)

					# add child to the population
					children[k] = bfs
	# if no crossover was performed, make a copy of the parents
	else:
		children = [[], []]
		for ind1 in list(parent1["individual"]):
			children[0] = children[0] + ind1
		for ind2 in list(parent2["individual"]):
			children[1] = children[1] + ind2

	return children


# Picks a node to add
def get_node_to_add(
		av_neighbors,
		loaded_data):

	deg = loaded_data["deg"]
	multiplex = loaded_data["multiplex"]

	# get the list of neighbors DEG
	de_neighbors = []
	for node in av_neighbors:
		if node in list(deg["gene"]):
			de_neighbors.append(node)

	# if there is at least one DE gene in the list of neighbors, keep
	if len(de_neighbors) > 0:
		av_neighbors = de_neighbors

	# calculate the number of times each node is neighbor
	incidences_ = Counter(av_neighbors)

	# get the nodes with highest incidence and keep nodes with a max incidence
	sorted_counts = incidences_.most_common()
	av_neighbors = []
	i = 0
	while sorted_counts[i][1] == sorted_counts[0][1] and i < len(sorted_counts):
		av_neighbors.append(sorted_counts[i][0])
		i += 1

	# pick a random node to add
	random_node = random.choice(av_neighbors)

	# get ID of the node to add
	new_node_id = get_id_of_nodes(random_node, multiplex[0])[0]

	return new_node_id


# Mutates a given list of nodes
def mutate_nodes(
		individual,
		ind_to_mut_net,
		nodes_to_mutate,
		pot_nodes_to_mutate,
		loaded_data):
	
	# get DEG
	deg = loaded_data["deg"]
	ind = []

	# if at least one of the nodes will be mutated
	# loop through the nodes to remove
	for i in nodes_to_mutate:
		mutated_network = ind_to_mut_net.delete_vertices(pot_nodes_to_mutate[i])

		# obtain neighbors of nodes
		neighbors_or_ind_ = []
		for node in ind_to_mut_net.vs:
			neighbors_or_ind_ = neighbors_or_ind_ + loaded_data["merged"].vs[loaded_data["merged"].neighbors(node)]["name"]

		neighbors_mut_ind_ = []
		for node in mutated_network.vs:
			neighbors_mut_ind_ = neighbors_mut_ind_ + loaded_data["merged"].vs[loaded_data["merged"].neighbors(node)]["name"]

		# delete nodes that originally belonged to the individual
		neighbors_or_ind = []
		for node in neighbors_or_ind_:
			if node not in ind_to_mut_net.vs["name"]:
				neighbors_or_ind.append(node)

		neighbors_mut_ind = []
		for node in neighbors_mut_ind_:
			if node not in ind_to_mut_net.vs["name"]:
				neighbors_mut_ind.append(node)

		# check if network is connected and there are available neighbors to add
		if mutated_network.is_connected() and len(neighbors_mut_ind) >0:
			# save changes
			ind = get_id_of_nodes(mutated_network.vs["name"], loaded_data["multiplex"][0])
			ind_to_mut_net = mutated_network
			av_neighbors = neighbors_mut_ind
		else:
			av_neighbors = neighbors_or_ind

		# if there is at least 1 available neighbor to be added
		if len(av_neighbors) > 0 and (ind and len(ind) < loaded_data["max_size"]):
			# pick a node to add
			new_node_id = get_node_to_add(av_neighbors, loaded_data)

			# add node to list and to network
			ind.append(new_node_id)
			ind_to_mut_net = loaded_data["merged"].induced_subgraph(ind)
		else:
			print("Attempt to add a new neighbor FAILED. Rolling back mutation")

	return ind


# Adds a node to a subnetwork
def add_node(
		ind,
		ind_to_mut_net,
		loaded_data):

	# get deg
	deg = loaded_data["deg"]

	# obtain neighbors of nodes
	av_neighbors_ = []
	for node in ind_to_mut_net.vs:
		av_neighbors_ = av_neighbors_ + loaded_data["merged"].vs[loaded_data["merged"].neighbors(node)]["name"]

	# delete from the list all the nodes that originally belonged to the ind
	av_neighbors = []
	for node in av_neighbors_:
		if node not in ind_to_mut_net.vs["name"]:
			av_neighbors.append(node)

	# if there is at least one available neighbor to be added
	if len(av_neighbors) > 0:
		# pick a node
		new_node_id = get_node_to_add(av_neighbors, loaded_data)

		# add node
		ind.append(new_node_id)
	else:
		print("No nodes to be mutated. Attempt to add a new neighbor FAILED")

	return ind


# Performs mutation
def mutation(
		individuals,
		multiplex,
		loaded_data):
	
	merged = loaded_data["merged"]
	deg = loaded_data["deg"]
	mutation_rate = loaded_data["mutation_rate"]

	mutants = [[]]*len(individuals)

	# loop through all the individuals to be mutated
	for i in range(len(individuals)):
		# generate a random number between 0 and 1
		p = np.random.uniform(0, 1)

		# check if mutation is to be performed
		if p <= mutation_rate:
			# make subnetwork
			ind_to_mut_net = merged.induced_subgraph(individual[i])
			# get nodes' degrees
			ind_to_mut_nodes_ = ind_to_mut_net.vs["name"]
			ind_to_mut_deg_ = ind_to_mut_net.degree(ind_to_mut_nodes)

			# remove DEG from the list
			ind_to_mut_nodes = []
			ind_to_mut_deg = []
			for j in range(len(ind_to_mut_deg_)):
				if ind_to_mut_nodes_[j] not in list(deg["gene"]):
					ind_to_mut_nodes.append(ind_to_mut_nodes_[j])
					ind_to_mut_deg.append(ind_to_mut_deg_[j])

			# get the list of nodes with the min degree (peripheral nodes) that can be mutated
			pot_nodes_to_mutate = []
			for j in range(len(ind_to_mut_nodes)):
				if ind_to_mut_deg[j] == min(ind_to_mut_deg):
					pot_nodes_to_mutate.append(ind_to_mut_nodes[j])

			# make vector with the nodes to mutate
			nodes_to_mutate_ = []
			for j in range(len(pot_nodes_to_mutate)):
				if np.random.uniform(0, 1) <= mutation_rate:
					nodes_to_mutate_.append(1)
				else:
					nodes_to_mutate_.append(0)

			# gets the list of chromosomes to be mutated
			nodes_to_mutate = []
			for j in range(len(nodes_to_mutate_)):
				if nodes_to_mutate_[j] == 1:
					nodes_to_mutate.append(nodes_to_mutate_.index(nodes_to_mutate_[j]))

			# if at least one of the nodes will be mutated
			if len(nodes_to_mutate) > 0:
				individuals[i] = mutate_nodes(individuals[i], ind_to_mut_net, nodes_to_mutate, pot_nodes_to_mutate, loaded_data)
			# if no nodes were removed, add node if max size allows it
			else:
				if len(individuals[i]) < loaded_data["max_size"]:
					individuals[i] = add_node(individuals[i], ind_to_mut_net, loaded_data)

		# save the individual in the mutants' list
		mutants[i] = individuals[i]

	return mutants


# 
def get_duplicated_inds(
		div_pop,
		threshold):
	
	# create all the combinations of 2 numbers with the ids of individuals
	index1 = []
	index2 = []
	js = []
	for i in range(len(div_pop)):
		for j in range(i+1, len(div_pop)):
			index1.append(i)
			index2.append(j)
			
			ind1 = div_pop.at[i, "individual"] 
			ind2 = div_pop.at[j, "individual"]

			# calculate the Jaccard similarity index = (intersection of A and B) / (union of A and B)
			js.append((len(list(set(ind1) & set(ind2))) / len(list(set(ind1) | set(ind2)))) * 100)

	sim = pd.DataFrame(list(itertools.zip_longest(index1, index2, js)), columns = ["index1", "index2", "js"]).fillna(0)

	# keep "duplicated" individuals
	sim = sim[sim["js"] >= threshold]

	return sim


# Replaces all duplicated individuals and those above the permitted threshold of similarity with random individuals
def replace_duplicated_inds(
		combined_population,
		loaded_data):
	
	div_pop = combined_population
	threshold = loaded_data["jaccard_similarity_threshold"]
	multiplex = loaded_data["multiplex"]
	# get duplicated individuals
	sim = get_duplicated_inds(div_pop, threshold)

	# verifiy if there are duplicated individuals
	if sim.shape[0] > 0:
		inds_to_remove = []

		# non-dominated sorting and crowding distance calculus
		sorted_pop = non_dom_sort(div_pop, loaded_data)

		i = 0
		# loop through all the duplicated individuals
		while i < sim.shape[0]:
			# ID of individuals
			ind1_id = sim.at[i, "index1"]
			ind2_id = sim.at[i, "index2"]

			if sim.at[i, "js"] < 100 and sorted_pop.at[ind1_id, "rank"] == 1 and sorted_pop.at[ind2_id, "rank"] == 1 and np.isinf(sorted_pop.at[ind1_id, "crowding_distance"]) and np.isinf(sorted_pop.at[ind2_id, "crowding_distance"]):
				print("Keeping similar individuals. Inf crowding distance.")
				print(sorted_pop.iloc[[ind1_id, ind2_id]])

			# tournament between the two individuals
			else:
				ind_to_keep = tournament_selection(2, sorted_pop.iloc[[ind1_id, ind2_id]])
				if ind_to_keep == sorted_pop.iloc[[ind1_id]]:
					inds_to_remove.append(ind2_id)
				else:
					inds_to_remove.append(ind1_id)

				# get and remove all future incidences of the ind
				ref_ = []
				for ind in inds_to_remove:
					if ind in list(sim["index1"]) and list(sim["index1"]).index(ind) not in ref_:
						ref_.append(list(sim["index1"]).index(ind))
					elif ind in list(sim["index2"]) and list(sim["index2"]).index(ind) not in ref_:
						ref_.append(list(sim["index2"]).index(ind))

				ref = []
				for r in ref_:
					if ref > i:
						ref.append(r)

				sim = sim.drop(ref)

			i += 1

		# remove the corresponding individuals
		div_pop = sorted_pop.drop(inds_to_remove)

	# generate as many new individuals as duplicated ones
	new_inds = generate_initial_pop(combined_population.shape[0] - div_pop.shape[0], multiplex, loaded_data)
	fitness_data = evaluate_population(new_inds, multiplex, loaded_data)

	div_pop = pd.DataFrame(data = {}, columns = ["individual"], dtype = object)

	for i in range(len(new_inds)):
		div_pop.at[i, "individual"] = new_inds[i]

	div_pop = pd.concat([div_pop, fitness_data], axis = 1)

	for i in range(len(new_inds)):
		div_pop.at[i, "rank"] = 0
		div_pop.at[i, "crowding_distance"] = 0

	div_pop = non_dom_sort(div_pop, loaded_data)

	return div_pop


# Performs replacement
def replacement(
		parents,
		children,
		loaded_data):
	
	pop_size = loaded_data["pop_size"]

	# combine the new population (offspring) with old population (parents)
	combined_population = pd.concat([parents, children])

	# replace duplicated individuals with random ones
	combined_population = replace_duplicated_inds(combined_population, loaded_data)

	# order combined population by rank
	combined_population = combined_population.sort_values(by = 'rank', ignore_index = True)

	# get last rank which will be in the replacement population
	last_rank = list(combined_population["rank"])[-1]

	# get last rank individuals
	last_rank_inds = combined_population[combined_population["rank"] == last_rank]

	# order by crowding distance
	last_rank_inds = last_rank_inds.sort_values(by = 'crowding_distance', ascending = False, ignore_index = True)

	# select new population for replacement
	new_population_for_replacement = combined_population[combined_population["rank"] < last_rank]

	i = 0
	while new_population_for_replacement.shape[0] < pop_size and i < last_rank_inds.shape[0]:
		new_population_for_replacement.append(last_rank_inds.iloc[[i]], ignore_index = True)
		i += 1

	return new_population_for_replacement


# Generate a new population from an existing one
def make_new_population(
		loaded_data,
		population):

	pop_size = loaded_data["pop_size"]
	max_number_of_attempts = loaded_data["max_number_of_attemps"]
	tournament_size = loaded_data["tournament_size"]
	multiplex = loaded_data["multiplex"]
	my_new_population = [[]]*pop_size # initialize vector

	# loop to generate the new population. In each loop, 2 children are created
	for i in range(0, pop_size, 2):
		attempts = 0 # counter for max. attemps to find parents
		keep_looking = True # flag to keep looking for parents

		while attempts < max_number_of_attempts and keep_looking:
			parent1 = tournament_selection(tournament_size, population) # selection
			parent2 = get_parent2(parent1, population, loaded_data)
			
			# verify if the search was unsuccessful
			if not parent2:
				# increment no. of attempts
				attempts = attempts + 1
			else:
				# parent 2 found - stop the search
				keep_looking = False

		# if unsuccessful search
		if attempts == max_number_of_attempts:
			print("Max. no. of attemps to find compatible parents")

			#generate two random individuals
			children = generate_initial_pop(2, multiplex, loaded_data)
		else:
			# mate parents
			children = crossover(parent1, parent2, loaded_data)
		
		# mutate children
		children = mutation(children, multiplex, loaded_data)
		# add child 1 to the population
		my_new_population[i] = children[0]
		# add child 2 to the population
		my_new_population[i+1] = children[1]

	# evaluate offspring
	fitness_data = evaluate_population(my_new_population, multiplex, loaded_data)

	# generate data frame with the individuals and their fitness
	new_population = pd.DataFrame(data = {}, columns = ["individual"], dtype = object)

	for i in range(len(my_new_population)):
		new_population.at[i, "individual"] = my_new_population[i]

	new_population = pd.concat([new_population, fitness_data], axis = 1)

	for i in range(len(my_new_population)):
		new_population.at[i, "rank"] = 0
		new_population.at[i, "crowding_distance"] = 0

	# prepare new population
	new_population_for_replacement = replacement(population, new_population, loaded_data)

	return new_population_for_replacement


# saves in a file the best pop_size individuals of the final population
def save_final_pop(
		best_inds_file,
		population,
		pop_size,
		network):
	
	# loop through the pop_size best individuals in the population
	for i in range(len(pop_size)):
		# get the individual's code
		ind = population.at[i, "individual"]

		# get the names of the corresponding nodes
		decoded_ind = network.vs["name"][ind]

		# save in files
		fd = open(best_inds_file, "a")
		fd.write(" ".join(decoded_ind) + "," + str(population.at[i, "average_nodes_score"]) + "," + str(population.at[i, "density"]) + "," + str(population.at[i, "rank"]) + "," + str(population.at[i, "crowding_distance"]))
		fd.close()


# Defines the function of the body of MOGAMUN
def mogamun_body(
		run_number,
		loaded_data,
		best_inds_path):

	best_inds_file = best_inds_path + "_Run_" + str(run_number) + ".txt"
	my_init_pop = generate_initial_pop(loaded_data["pop_size"], loaded_data["multiplex"], loaded_data)
	fitness_data = evaluate_population(my_init_pop, loaded_data["multiplex"], loaded_data)
	population = pd.DataFrame(data = {}, columns = ["individual"], dtype = object)

	for i in range(len(my_init_pop)):
		population.at[i, "individual"] = my_init_pop[i]

	population = pd.concat([population, fitness_data], axis = 1)

	# obtain ranking and crowding distances
	population = non_dom_sort(population, loaded_data)

	# initilizes the number of generation
	g <- 1
	
	generation = []
	best_average_nodes_score = []
	best_density = []

	if_all_rank = True
	for rank in list(population["rank"]):
		if rank != 1:
			if_all_rank = False
			break

	# evolution's loop for g generations or until all inds have rank = 1
	while g <= loaded_data["generations"] and not if_all_rank:
		population = make_new_population(loaded_data, population)

		# add the best values for the two objective functions
		generation.append(g)
		best_average_nodes_score.append(max(list(population["average_nodes_score"])))
		best_density.append(max(list(population["density"])))

		print("Run " + str(run_number) + ". Gen. " + str(g) + " completed")

		# increments the generation
		g += 1

	# saves data in files
	stats_gen = pd.DataFrame(list(itertools.zip_longest(generation, best_average_nodes_score, best_density)), columns = ["generation", "best_average_nodes_score", "best_density"]).fillna(0)

	stats_gen.to_csv(path_or_buf = best_inds_path + "StatisticsPerGeneration_Run" + str(run_number) + ".csv", index = False)

	save_final_pop(best_inds_file, population, loaded_data["pop_size"], loaded_data["multiplex"][0])

	print("FINISH TIME, RUN " + str(run_number) + ": " + datetime.now().strftime("%H:%M:%S"))