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
    nodes_scores
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