import argparse
import json
import os

import numpy as np
from models.sup_models import load_data, split_data, ClassicalModel, NeuralNetworkModel, preprocess_data
import logging
from datetime import datetime
# Import necessary modules or functions for model handling and data processing
# from your_project_module import load_data, split_data, preprocess_data, ClassicalModel, NeuralNetworkModel

def main(config, model_filepath=None, action="train"):
    """
    Main function to handle training, saving, loading, and evaluating models based on the specified action.

    Args:
        config (dict): Configuration dictionary for data handling, model initialization, and training.
        model_filepath (str, optional): File path to save or load the model.
        action (str): Action to perform - 'train', 'save', 'load', or 'evaluate'.

    Returns:
        None
    """
    # Load Data
    data_path = "C:\\Users\\gqu\\OneDrive - UTHealth Houston\\projects\\Genevic\\data"
    # data_path = r"/mnt/c/Users/Gang/OneDrive - UTHealth Houston/projects/Genevic/data"
    data_path = r"/mnt/c/Users/gqu/Downloads/data"
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data directory not found at path: {data_path}")

    X, y, rsID = load_data(data_path)
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Create the plot using seaborn
    plt.figure(figsize=(8, 6))
    sns.histplot(y, kde=True, color='skyblue', edgecolor='black', bins=30)

    # Add titles and labels
    plt.title("Distribution of y")
    plt.xlabel("y values")
    plt.ylabel("Frequency")

    # Display the plot
    plt.show()

    # keep the top N
    n = 20

    # Get the original indices (0-based) of the data
    indices = list(range(len(y)))

    # Sort indices based on the y values (in descending order)
    sorted_indices = sorted(indices, key=lambda i: y[i], reverse=True)

    # Get the indices of the top n samples
    top_n_indices = set(sorted_indices[:n])

    # Extract top n samples (sorted by descending y)
    top_n_X = [X[i] for i in sorted_indices[:n]]
    top_n_y = [y[i] for i in sorted_indices[:n]]
    top_n_rsID = [rsID[i] for i in sorted_indices[:n]]

    # Extract the remaining samples in their original order
    remaining_X = np.array([X[i] for i in indices if i not in top_n_indices])
    remaining_y = np.array([y[i] for i in indices if i not in top_n_indices]).reshape(-1,)
    remaining_rsID = [rsID[i] for i in indices if i not in top_n_indices]




    # logging.info(top_n_rsID)
    # print("Data successfully loaded.")
    logging.info("Data successfully loaded.")
    # Split Data
    test_size = config.get('test_size', 0.2)
    random_state = config.get('random_state', 42)
    X_train, X_test, y_train, y_test, rsID_train, rsID_test = split_data(
        remaining_X, remaining_y, remaining_rsID, test_size=test_size
    )

    X_train_scaled, X_test_scaled, y_train_transformed, y_test_transformed, X_scaled, y_transformed = preprocess_data(
        X_train, X_test, y_train, y_test, remaining_X, remaining_y)
    # print(f"Data split completed: {len(X_train)} training samples, {len(X_test)} testing samples.")
    logging.info("Data split completed: {len(X_train)} training samples, {len(X_test)} testing samples.")
    # Model Selection
    method = config.get('method', 'random_forest')
    model_class = ClassicalModel if method in ['linear_regression', 'decision_tree', 'random_forest', 'svr', 'lasso',
                                               'ridge'] else NeuralNetworkModel
    model = model_class(config)

    # Perform actions based on user input
    if action == "train":
        # model.train(X_train, y_train)
        model.train(X_train_scaled, y_train_transformed)
        # print("Training completed successfully.")
        logging.info("Training completed successfully.")
        # Save the trained model if a filepath is provided
        if model_filepath:
            model.save_model(model_filepath)
            print(f"Model saved to {model_filepath}")
        # results = model.evaluate(X_test, y_test)
        results = model.evaluate(X_test_scaled, y_test_transformed)
        # print(results)
        logging.info(results)
        predictions, ranks, rsID_to_rank, ordered_rsIDs = model.predict(X_scaled, y_transformed, rsID)
        logging.info(top_n_rsID+ordered_rsIDs)
    elif action == "load" and model_filepath:
        model = model_class.load_model(model_filepath, config)
        # print(f"Model loaded from {model_filepath}")
        logging.info("Model loaded from {model_filepath}")
        # Evaluate the loaded model
        results = model.evaluate(X_test, y_test)
        # print(f"Evaluation Results: {results}")
        logging.info("Evaluation Results: {results}")
    elif action == "evaluate":
        results = model.evaluate(X_test, y_test)
        # print(f"Evaluation Results: {results}")
        logging.info("Evaluation Results: {results}")

    else:
        raise ValueError("Invalid action specified. Choose from 'train', 'save', 'load', or 'evaluate'.")


def parse_arguments():
    """
    Parse command-line arguments for the main script.

    Returns:
        Namespace: Parsed arguments including config file path, model file path, and action.
    """
    parser = argparse.ArgumentParser(description="Train, Save, Load, and Evaluate Machine Learning Model")
    parser.add_argument('--config', type=str, required=True, help="Path to JSON configuration file.")
    parser.add_argument('--model_filepath', type=str, default=None, help="Path to save or load the model file.")
    parser.add_argument('--action', type=str, choices=['train', 'save', 'load', 'evaluate'], default="train",
                        help="Action to perform: 'train', 'save', 'load', or 'evaluate'.")
    return parser.parse_args()


if __name__ == "__main__":
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    if not os.path.exists("logs"):
        os.makedirs("logs")
    # Build your time-based log filename
    logging.basicConfig(
        level=logging.DEBUG,
        filename=f"logs/logfile_{current_time}.log",  # Update to your desired path
        filemode="a+",  # 'a+' appends logs; use 'w' to overwrite each run
        format="%(asctime)-15s %(levelname)-8s %(message)s"
    )
    args = parse_arguments()
    # Load configuration from the specified JSON file
    if not os.path.isfile(args.config):
        raise FileNotFoundError(f"Configuration file not found at path: {args.config}")

    with open(args.config, 'r') as f:
        config = json.load(f)

    # print(f"Configuration loaded from {args.config}.")
    logging.info("Configuration loaded from {args.config}.")
    logging.info(config)

    try:
        main(config, model_filepath=args.model_filepath, action=args.action)
    except Exception as e:
        # print(f"An error occurred during execution: {e}")
        logging.info("An error occurred during execution: {e}")
