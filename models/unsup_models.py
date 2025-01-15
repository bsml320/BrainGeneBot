import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import kendalltau, rankdata
import random


def flatten_ranks(ranks):
    """
    Flatten ranks if they are lists.
    """
    flattened_ranks = []
    for rank in ranks:
        if isinstance(rank, list):
            flattened_ranks.extend(rank)
        else:
            flattened_ranks.append(rank)
    return flattened_ranks


def handle_rank_value(rank):
    """
    Handle rank values that could be lists by averaging the values.
    """
    if isinstance(rank, list):
        return np.mean(rank)
    return rank


def borda_method(rankings, method):
    """
    Borda method for rank aggregation.
    Method: MEAN, MED, GEO, L2
    - MEAN: Takes the average of the ranks.
    - MED: Takes the median of the ranks.
    - GEO: Uses the geometric mean of the ranks.
    - L2: Uses the L2 norm to aggregate ranks.
    """
    methods = {'mean': np.mean, 'med': np.median, 'geo': lambda x: np.prod(x) ** (1 / len(x)),
               'l2': lambda x: np.sqrt(np.sum(np.square(x)) / len(x))}
    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
    aggregated_scores = defaultdict(list)

    for ranking in rankings:
        for item in all_items:
            if item in ranking:
                rank_value = handle_rank_value(ranking[item])
                aggregated_scores[item].append(rank_value)
            else:
                aggregated_scores[item].append(np.nan)

    aggregated_rank = {item: methods[method.lower()](pd.Series(scores).dropna()) for item, scores in
                       aggregated_scores.items()}
    # Convert aggregated ranks into scores using 'dense' ranking method
    scores_list = [aggregated_rank[item] for item in all_items]
    rank_scores = rankdata(scores_list, method='dense')

    # Create a pandas Series from the rank scores, with items as index
    return pd.Series(rank_scores, index=all_items)


def markov_chain(rankings, method):
    """
    Markov Chain methods for rank aggregation.
    Method: MC1, MC2, MC3, MC4
    - MC1: Basic Markov Chain model using pairwise comparisons.
    - MC2: Weighted Markov Chain to give certain rankings more influence.
    - MC3: Rank-dependent Markov Chain, where transition probabilities are adjusted based on rank differences.
    - MC4: Random walk with restart, similar to the PageRank algorithm, to prevent getting stuck in cycles.
    """
    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
    num_items = len(all_items)
    if num_items == 0:
        raise ValueError("No items to rank.")

    item_index = {item: idx for idx, item in enumerate(all_items)}
    transition_matrix = np.zeros((num_items, num_items))

    for ranking in rankings:
        sorted_items = sorted(ranking.keys(), key=lambda x: ranking[x])
        for i in range(len(sorted_items)):
            for j in range(i + 1, len(sorted_items)):
                transition_matrix[item_index[sorted_items[i]], item_index[sorted_items[j]]] += 1

    # Normalizing transition matrix
    transition_matrix += np.eye(num_items)  # Avoid division by zero
    transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)

    if method == 'mc2':
        transition_matrix *= 1.2
        transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)
    elif method == 'mc3':
        for ranking in rankings:
            sorted_items = sorted(ranking.keys(), key=lambda x: ranking[x])
            for i in range(len(sorted_items)):
                for j in range(i + 1, len(sorted_items)):
                    rank_diff = abs(ranking[sorted_items[j]] - ranking[sorted_items[i]])
                    transition_matrix[item_index[sorted_items[i]], item_index[sorted_items[j]]] += rank_diff
        transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)
    elif method == 'mc4':
        restart_prob = 0.15
        transition_matrix = restart_prob / num_items + (1 - restart_prob) * transition_matrix

    stationary_distribution = np.ones(num_items) / num_items
    prev_distribution = np.zeros(num_items)
    for _ in range(100):
        if np.allclose(stationary_distribution, prev_distribution, atol=1e-6):
            break
        prev_distribution = stationary_distribution
        stationary_distribution = stationary_distribution @ transition_matrix

    rank_scores = rankdata(stationary_distribution, method='dense')
    return pd.Series(rank_scores, index=all_items)


def thurstone_model(rankings):
    """
    Thurstone's model for rank aggregation.
    Assumes higher score values should translate to lower preference ranks.
    """
    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
    scores = defaultdict(float)
    count_scores = defaultdict(int)

    for ranking in rankings:
        # Use rank values as they are to compute average scores
        ranked_items = [ranking[item] for item in ranking.keys()]
        for item, score in zip(ranking.keys(), ranked_items):
            scores[item] += score
            count_scores[item] += 1

    # Average the scores for each item
    for item in scores:
        if count_scores[item] > 0:
            scores[item] /= count_scores[item]

    # Convert average scores to ranks where a higher score gives a higher rank number
    sorted_items_by_score = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    ranks = rankdata([score for item, score in sorted_items_by_score], method='dense')

    # Create a pandas Series from the ranks
    return pd.Series(ranks, index=[item for item, score in sorted_items_by_score])


def stuart_method(rankings):
    """
    Stuart method for rank aggregation.
    Aggregates ranks by summing the original rank values from all rankings.
    Uses rankdata to assign dense ranks to the aggregated scores.
    Lower sum of ranks indicates better performance.
    """
    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
    aggregated_rank = defaultdict(int)

    for ranking in rankings:
        for item, rank in ranking.items():
            aggregated_rank[item] += rank  # Sum the rank values directly

    # Applying rankdata to assign ranks based on aggregated ranks
    # Lower aggregated ranks get lower (better) rank numbers
    rank_scores = rankdata([aggregated_rank[item] for item in all_items], method='dense')
    return pd.Series(rank_scores, index=all_items)

def rra_method(rankings):
    """
    Robust Rank Aggregation (RRA) method.
    RRA aims to assess the significance of items being ranked higher across multiple rankings.
    """
    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
    num_items = len(all_items)
    item_index = {item: idx for idx, item in enumerate(all_items)}
    p_values = np.ones(num_items)

    for ranking in rankings:
        ranked_items = rankdata([handle_rank_value(ranking[item]) for item in ranking.keys()])
        for item, score in zip(ranking.keys(), ranked_items):
            p_values[item_index[item]] *= score / num_items

    return pd.Series(p_values, index=all_items)


def bayesian_method(rankings, method):
    """
    Bayesian method for rank aggregation.
    Method: BARD, BIRRA, BIG
    - BARD: Bayesian rank aggregation with negative scores.
    - BIRRA: Bayesian rank aggregation with squared scores.
    - BIG: Bayesian rank aggregation with logarithmic scores.
    """
    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
    scores = defaultdict(float)

    # Aggregate scores from rank data using ordinal ranks directly from the values provided
    for ranking in rankings:
        ranked_items = rankdata([ranking[item] for item in ranking.keys()], method='ordinal')
        for item, score in zip(ranking.keys(), ranked_items):
            scores[item] += score

    # Transform scores according to the selected method
    if method == 'bard':
        transformed_scores = {item: scores[item] for item in all_items}
    elif method == 'birra':
        transformed_scores = {item: scores[item]**2 for item in all_items}
    elif method == 'big':
        transformed_scores = {item: np.log1p(scores[item]) for item in all_items}
    else:
        raise ValueError(f"Bayesian method {method} is not supported.")

    # Convert transformed scores into a list matching the sorted order of all_items
    scores_list = [transformed_scores[item] for item in all_items]

    # Apply rankdata with 'dense' to assign ranks based on transformed scores
    rank_scores = rankdata(scores_list, method='dense')

    # Create a pandas Series from the rank scores, with items as index
    return pd.Series(rank_scores, index=all_items)


def stochastic_optimization(rankings, method, seed=42, iterations=10000):
    """
    Stochastic optimization for rank aggregation with reproducibility.
    Method: CEMC (Kendall), CEMC (Spearman's footrule)
    - CEMC (Kendall): Minimizes Kendall's tau distance between rankings.
    - CEMC (Spearman's footrule): Minimizes Spearman's footrule distance between rankings.
    """
    random.seed(seed)  # Seed the random number generator for reproducibility
    np.random.seed(seed)

#    all_items = sorted(set(item for ranking in rankings for item in ranking.keys()))
#     all_items = sorted(set(item for ranking in rankings for item in ranking),
#                        key=lambda x: np.mean([ranking.get(x, float('inf')) for ranking in rankings]))
    all_ranks = [rank for ranking in rankings for rank in ranking.values()]
    all_items = list(set([item for ranking in rankings for item in ranking.keys()]))
    num_items = len(all_items)
    best_rank = np.arange(1, num_items + 1)
    best_score = np.inf

    for _ in range(iterations):
        candidate_rank = np.random.permutation(best_rank)
        current_score = 0

        if method == 'cemc_kendall':
            for ranking in rankings:
                rank_vector = np.array([ranking.get(item, num_items + 1) for item in all_items])
                tau, _ = kendalltau(candidate_rank, rank_vector)
                current_score -= tau
        elif method == 'cemc_spearman':
            for ranking in rankings:
                rank_vector = np.array([ranking.get(item, num_items + 1) for item in all_items])
                current_score += np.sum((candidate_rank - rank_vector) ** 2)
        else:
            raise ValueError(f"Stochastic optimization method {method} is not supported.")

        if current_score < best_score:
            best_rank, best_score = candidate_rank, current_score

    rank_scores = rankdata(best_rank, method='dense')
    return pd.Series(rank_scores, index=all_items)

def rank_aggregation(data, config):
    """
    Aggregate ranks from multiple studies using various rank aggregation methods.
    Args:
    - data: dict of study_id: pd.DataFrame with columns [gene_variant_id, effect_weight, ranks]
    - config: dict of configuration settings
    Returns: pd.DataFrame with columns [gene_variant_id, aggregated_rank]
    """
    rankings = []
    for study_id, study_df in data.items():
        study_df['hm_rsID'] = study_df['hm_rsID'].astype(str)
        rankings.append(dict(zip(study_df['hm_rsID'], [handle_rank_value(rank) for rank in study_df['ranks']])))

    method = config.get('method').lower()

    print(method)
    if method in ['mean', 'med', 'geo', 'l2']:
        aggregated_rank = borda_method(rankings, method)

    elif method in ['mc1', 'mc2', 'mc3', 'mc4']:
        aggregated_rank = markov_chain(rankings, method)
    elif method == 'thurstone':
        aggregated_rank = thurstone_model(rankings)
    elif method == 'stuart':
        aggregated_rank = stuart_method(rankings)
    elif method == 'rra':
        aggregated_rank = rra_method(rankings)
    elif method in ['bard', 'birra', 'big']:
        aggregated_rank = bayesian_method(rankings, method)
    elif method in ['cemc_kendall', 'cemc_spearman']:
        aggregated_rank = stochastic_optimization(rankings, method)
    else:
        raise ValueError(f"Method {method} is not supported.")

    aggregated_rank.index = aggregated_rank.index.map(str)
    result_df = pd.DataFrame({'gene_variant_id': aggregated_rank.index, 'aggregated_rank': aggregated_rank.values})
    return result_df


def process_and_aggregate(data_configs):
    """
    Processes multiple data sets with specified methods and aggregates into a single DataFrame.

    :param data_configs: List of tuples, each containing a DataFrame and its config dict.
    :return: A single DataFrame with aggregated ranks across specified methods.
    """
    processed_data_list = []

    for data, config in data_configs:
        method_name = config['method']
        # Assuming 'data' is a DataFrame with at least 'gene_variant_id' and 'score'
        processed_data = data[['gene_variant_id', 'score']].rename(columns={'score': method_name})
        processed_data_list.append(processed_data)

    # Start with the first DataFrame
    aggregated_df = processed_data_list[0]

    # Iteratively merge each DataFrame into the aggregated DataFrame
    for additional_data in processed_data_list[1:]:
        aggregated_df = pd.merge(aggregated_df, additional_data, on='gene_variant_id', how='inner')

    return aggregated_df


# if __name__ == "__main__":
#     import pickle
#
#     # Example Usage:
#     data = {
#         'study1': pd.DataFrame(
#             {'hm_rsID': ['1', '2', '3'], 'effect_weight': [0.1, 0.2, 0.15], 'ranks': [1, 2, 3.5]}),
#         'study2': pd.DataFrame(
#             {'hm_rsID': ['1s', '2', '4'], 'effect_weight': [0.3, 0.25, 0.4], 'ranks': [[2, 3], 3, 1]}),
#     }
#
#     with open(r'C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\data_dict.pkl', 'rb') as f:
#         data = pickle.load(f)
#         print(data['PGS003953'])
#     # Iterate through each study/dataframe in the data dictionary
#     for study, df in data.items():
#         print(f"Checking {study}:")
#         # Check which columns contain NaNs
#         for column in df.columns:
#             if df[column].isnull().any():
#                 nan_count = df[column].isnull().sum()
#                 print(f"  -> {column} contains {nan_count} NaN(s)")
#     config = {'method': 'mean'}
#
#     aggregated_ranks_df = rank_aggregation(data, config)
#     print(aggregated_ranks_df)
