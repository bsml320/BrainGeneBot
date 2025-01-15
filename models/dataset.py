import pandas as pd
import numpy as np
from tqdm import tqdm

import numpy as np
import pandas as pd
from tqdm import tqdm


def feature_engineering(data_dict, labels_df):
    all_features = []
    all_labels = []
    all_variant_ids = []

    # Convert label data to a dictionary for fast lookup
    label_dict = labels_df.set_index('variant_id')['priority_score'].fillna(0).to_dict()

    # Combine all DataFrames into one DataFrame
    combined_df = pd.concat(data_dict.values(), ignore_index=True)

    # Process grouped data by hm_rsID
    for hm_rsID, group in tqdm(combined_df.groupby('hm_rsID'), desc="Processing hm_rsIDs"):
        effect_weights = group['effect_weight'].values


        # Calculate effect statistics if sufficient data is available
        if effect_weights.size > 0:
            mean_effect = np.mean(effect_weights)
            sum_effect = np.sum(effect_weights)
            min_effect = np.min(effect_weights)
            max_effect = np.max(effect_weights)
            range_effect = np.ptp(effect_weights)
            std_effect = np.std(effect_weights, ddof=1) if effect_weights.size > 1 else 0
        else:
            mean_effect = sum_effect = min_effect = max_effect = range_effect = std_effect = np.nan

        # Handling ranks: concatenate all rank arrays/lists
        ranks = [np.array(x) if isinstance(x, list) else np.array([x]) for x in group['ranks'] if not pd.isna(x)]

        # Prepare and calculate rank features if data is available
        if ranks:
            ranks_array = np.concatenate(ranks)
            rank_mean = np.mean(ranks_array) if ranks_array.size > 0 else np.nan
            rank_sum = np.sum(ranks_array) if ranks_array.size > 0 else np.nan
            rank_range = np.ptp(ranks_array) if ranks_array.size > 1 else 0
            rank_std = np.std(ranks_array, ddof=1) if ranks_array.size > 1 else 0
            rank_count = len(ranks_array)
        else:
            rank_mean = rank_sum = rank_range = rank_std = rank_count = np.nan

        features = np.array([mean_effect, sum_effect, min_effect, max_effect, range_effect, std_effect,
                             rank_mean, rank_sum, rank_range, rank_std, rank_count])
        all_features.append(features)
        all_labels.append(label_dict.get(hm_rsID, 0))
        all_variant_ids.append(hm_rsID)

    return np.array(all_features), np.array(all_labels), all_variant_ids




if __name__ == "__main__":
    import pickle
    # Load data
    labels_df = pd.read_csv(r"C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\AD_GWAS_Priority_Scores.csv")
    with open(r'C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\\data_dict.pkl', 'rb') as f:
        data_dict = pickle.load(f)

    # Create the dataset
    X, y, all_variant_ids = feature_engineering(data_dict, labels_df)
    #
    print("Features (X):")
    print(X, X.shape)
    print("\nLabels (y):")
    print(y, y.shape)
    print(all_variant_ids)
    print(len(set(all_variant_ids)))
    np.savez(r'C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data\\data.npz', allow_pickle=True)
