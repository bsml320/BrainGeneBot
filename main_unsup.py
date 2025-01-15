import argparse
import pickle
import pandas as pd
from dask.cli import config_get
from models.unsup_models import rank_aggregation
import logging
from datetime import datetime
import os


current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
def parse_arguments():
    """
    Parse command-line arguments for the data path and optionally an output file path.
    """
    parser = argparse.ArgumentParser(description="Perform rank aggregation.")
    parser.add_argument('--data_path', type=str, required=True,
                        help="Path to the pickle file containing the data.")

    parser.add_argument('--config', type=str, default='mean',
                        help="method")

    parser.add_argument('--prefix', type=str, default=str(current_time),
                        help="prefix")
    parser.add_argument('--output_path', type=str, default='logs/unsup/',
                        help="Path where the aggregated rank results will be saved. Defaults to 'aggregated_ranks.csv'.")
    return parser.parse_args()


def main(data_path, config, output_path, prefix):
    """
    Main function to perform rank aggregation and save the results to a file.

    Args:
        data_path (str): Path to the pickle file containing the data.
        output_path (str): Path to save the aggregated rank results.
    """
    # Load data from the specified pickle file
    with open(data_path, 'rb') as f:
        data = pickle.load(f)


    # Iterate through each study/dataframe in the data dictionary
    for study, df in data.items():
        logging.info(f"Checking {study}:")
        # Check which columns contain NaNs
        for column in df.columns:
            if df[column].isnull().any():
                nan_count = df[column].isnull().sum()
                print(f"  -> {column} contains {nan_count} NaN(s)")

    config_g = {'method': config}
    # Perform rank aggregation using the loaded data and the hard-coded configuration
    aggregated_ranks_df = rank_aggregation(data, config_g)
    print(aggregated_ranks_df)

    # Save the aggregated ranks to the specified output file
    aggregated_ranks_df.to_csv(os.path.join(output_path,prefix+'_'+'aggregated_ranks'+str(current_time)+'.csv'), index=False)
    logging.info(f"Aggregated ranks saved to {output_path}")


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.DEBUG,
        filename=f"logs/unsup/logfile_{current_time}.log",  # Update to your desired path
        filemode="a+",  # 'a+' appends logs; use 'w' to overwrite each run
        format="%(asctime)-15s %(levelname)-8s %(message)s"
    )

    args = parse_arguments()
    config = {'method':'mean'}
    # logging.info(config)
    main(args.data_path, args.config, args.output_path, args.prefix)
