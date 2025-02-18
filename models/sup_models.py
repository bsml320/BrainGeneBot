import numpy as np
import joblib
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler,  PowerTransformer
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression, Lasso, Ridge
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.compose import TransformedTargetRegressor
import pickle, os
import sys
from models.dataset import feature_engineering
import random
from Bio import Entrez
import matplotlib.pyplot as plt
import seaborn as sns


def rsid2gene(rsids):
    Entrez.email = 'adamrt9319@gmail.com'
    names = []
    for snp_id in rsids:
        numeric_id = snp_id.replace('rs', '')
        record = Entrez.read(Entrez.elink(dbfrom="snp",
                                          id=numeric_id,
                                          db="gene"))
        # Check if LinkSetDb and links are present
        if 'LinkSetDb' in record[0] and len(record[0]['LinkSetDb']) > 0:
            links = record[0]['LinkSetDb'][0]['Link']
            for gene_link in links:
                gene_id = gene_link['Id']
                handle = Entrez.esummary(db="gene", id=gene_id)
                uid_record = Entrez.read(handle)
                handle.close()

                uid_summary = uid_record["DocumentSummarySet"]['DocumentSummary'][0]
                names.append(uid_summary['Name'])
        else:
            print(f"No gene links found for {snp_id}")
    return names

def seed_it(seed):
    """
    set random seed for reproducibility
    :param seed: seed number
    :type seed: int
    :return:
    :rtype:
    """
    random.seed(seed)
    os.environ["PYTHONSEED"] = str(seed)
    np.random.seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.enabled = True
    torch.manual_seed(seed)

def load_data(data_path):
    # Assuming data is available as a CSV or in a DataFrame
    # Replace this part with the actual data loading process if necessary

    labels_df = pd.read_csv(os.path.join(data_path,"AD_GWAS_Priority_Scores.csv"))
    with open(os.path.join(data_path,"data_dict.pkl"), 'rb') as f:
        data_dict = pickle.load(f)

    X, y, all_variant_ids = feature_engineering(data_dict, labels_df)
    return X, y, all_variant_ids


def split_data(X, y, all_variant_ids, test_size=0.2):
    X_train, X_test,  y_train, y_test, rsID_train, rsID_test = train_test_split(X, y,all_variant_ids, test_size=test_size, shuffle=False)
    return X_train, X_test, y_train, y_test, rsID_train, rsID_test



def preprocess_data(X_train, X_test, y_train, y_test, X=None, y=None):
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # X_train_scaled = X_train
    # X_test_scaled = X_test
    power_transformer = PowerTransformer(method='yeo-johnson')
    power_transformer.fit(y_train.reshape(-1, 1))
    y_train_transformed = power_transformer.transform(y_train.reshape(-1, 1))
    y_test_transformed = power_transformer.transform(y_test.reshape(-1, 1))

    if X is not None and y is not None:
        X_scaled = scaler.transform(X)

        y_transformed = power_transformer.transform(y.reshape(-1, 1))

        return X_train_scaled, X_test_scaled, y_train_transformed, y_test_transformed, X_scaled, y_transformed
    return X_train_scaled, X_test_scaled, y_train_transformed, y_test_transformed


class BaseModel:
    def train(self, X_train, y_train, rsID_train=None):
        """
        Train the model on the provided training data.

        Args:
            X_train: Training feature data.
            y_train: Training labels.
            rsID_train: Optional identifier data.
        """
        raise NotImplementedError("Train method not implemented.")

    def evaluate(self, X_test, y_test, rsID_test=None):
        """
        Evaluate the model on test data and return performance metrics.

        Args:
            X_test: Test feature data.
            y_test: Test labels.
            rsID_test: Optional identifier data.
        """
        raise NotImplementedError("Evaluate method not implemented.")

    def save_model(self, filepath):
        """
        Save the model to a specified file path.

        Args:
            filepath (str): Path to the file where the model should be saved.
        """
        raise NotImplementedError("Save_model method not implemented.")

    def predict(self, X, y, rsID_test=None):
        """
        predict the results on data and return performance metrics.

        Args:
            X: Feature data.
            y: Test labels.
            rsID_test: Optional identifier data.
        """
        raise NotImplementedError("Save_model method not implemented.")

    @classmethod
    def load_model(cls, filepath, config=None):
        """
        Load a model from a file and return an instance of the class.

        Args:
            filepath (str): Path to the file from which the model should be loaded.
            config (dict, optional): Configuration dictionary for model initialization.

        Returns:
            An instance of the model class with the loaded state.
        """
        raise NotImplementedError("Load_model method not implemented.")

# Implement the Neural Network Model
class NeuralNetworkModel(BaseModel):
    def __init__(self, config):
        """
        Initialize the neural network model based on the configuration dictionary.

        Args:
            config (dict): Configuration dictionary containing model parameters.
        """
        self.input_dim = config.get('input_dim', 11)
        self.learning_rate = config.get('learning_rate', 0.001)
        self.epochs = config.get('epochs', 50)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model = self.define_model().to(self.device)
        self.criterion = nn.MSELoss()
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.learning_rate)

    def define_model(self):
        """
        Define the architecture of the neural network.

        Returns:
            nn.Sequential: Neural network model architecture.
        """
        return nn.Sequential(
            nn.Linear(self.input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
        )

    def train(self, X_train, y_train, rsID_train=None):
        """
        Train the model on the provided training data.

        Args:
            X_train: Training feature data.
            y_train: Training labels.
            rsID_train: Annotation list
        """
        X_train = torch.tensor(X_train, dtype=torch.float32).to(self.device)
        y_train = torch.tensor(y_train, dtype=torch.float32).view(-1, 1).to(self.device)
        self.model.train()

        for epoch in range(self.epochs):
            self.optimizer.zero_grad()
            outputs = self.model(X_train)
            loss = self.criterion(outputs, y_train)
            loss.backward()
            self.optimizer.step()

            if (epoch + 1) % 10 == 0:
                print(f'Epoch [{epoch + 1}/{self.epochs}], Loss: {loss.item():.4f}')

    def evaluate(self, X_test, y_test, rsID_test=None):
        """
        Evaluate the model on test data and return performance metrics.

        Args:
            X_test: Test feature data.
            y_test: Test labels.
            rsID_test: Annotation list

        Returns:
            dict: Dictionary containing Mean Squared Error and R2 Score.
        """
        X_test = torch.tensor(X_test, dtype=torch.float32).to(self.device)
        self.model.eval()
        with torch.no_grad():
            predictions = self.model(X_test).cpu().numpy()
        mse = mean_squared_error(y_test, predictions)
        r2 = r2_score(y_test, predictions)
        return {'Mean Squared Error': mse, 'R2 Score': r2}

    def predict(self, X, y, rsID=None):
        X = torch.tensor(X, dtype=torch.float32).to(self.device)
        self.model.eval()
        with torch.no_grad():
            predictions = self.model(X).cpu().numpy()

        sorted_indices = np.argsort(predictions)

        ranks = np.argsort(sorted_indices) + 1

        if rsID is not None:
            rsID = np.array(rsID)
            # Create a dictionary mapping each rsID to its corresponding rank
            rsID_to_rank = dict(zip(rsID, ranks))
            # Sort rsIDs based on their ranks (ascending order)
            ordered_rsIDs = [x for _, x in sorted(zip(ranks, rsID), key=lambda pair: pair[0])]
            return predictions, ranks, rsID_to_rank, ordered_rsIDs

        return predictions, ranks
    def save_model(self, filepath):
        """
        Save the trained model's state dictionary to a file.

        Args:
            filepath (str): Path to the file where the model should be saved.
        """
        torch.save(self.model.state_dict(), filepath)
        print(f"Model state dictionary saved to {filepath}")

    @classmethod
    def load_model(cls, filepath, config=None):
        """
        Load a model's state dictionary from a file and return an instance of NeuralNetworkModel.

        Args:
            filepath (str): Path to the file from which the model state dictionary should be loaded.
            config (dict): Configuration dictionary for initializing the model structure.

        Returns:
            NeuralNetworkModel: An instance of NeuralNetworkModel with the loaded state dictionary.
        """
        """
        Load a model's state dictionary from a file and return an instance of NeuralNetworkModel.

        Args:
            filepath (str): Path to the file from which the model state dictionary should be loaded.
            config (dict, optional): Configuration dictionary for initializing the model structure.

        Returns:
            NeuralNetworkModel: An instance of NeuralNetworkModel with the loaded state dictionary.
        """
        if config is None:
            raise ValueError("Configuration dictionary must be provided to initialize the model structure.")

        # Initialize model instance with configuration
        instance = cls(config)
        # Load the state dictionary
        instance.model.load_state_dict(torch.load(filepath, map_location=instance.device))
        instance.model.to(instance.device)
        print(f"Model loaded from {filepath}")
        return instance

# Implement Classical Machine Learning Models with Hyperparameter Configuration
class ClassicalModel(BaseModel):
    def __init__(self, config):
        """
        Initialize the model based on the configuration dictionary.

        Args:
            config (dict): Configuration dictionary containing model type and parameters.
        """
        method = config.get('method', 'linear_regression')
        model_params = config.get('model_params', {})
        use_ttr = config.get('use_transformed_target_regressor', False)
        models = {
            'linear_regression': LinearRegression,
            'decision_tree': DecisionTreeRegressor,
            'random_forest': RandomForestRegressor,
            'svr': SVR,
            'lasso': Lasso,
            'ridge': Ridge,
            # ... add other models
        }
        ModelClass = models.get(method)
        if ModelClass is None:
            raise ValueError(f"Unsupported method: {method}")
        base_model = ModelClass(**model_params)
        if use_ttr:
            # Apply exponential transformation to ensure positive predictions
            self.model = TransformedTargetRegressor(
                regressor=base_model,
                func=np.log1p,   # Apply log(1 + y)
                inverse_func=np.expm1  # Apply exp(y) - 1
            )
        else:
            self.model = base_model

    def train(self, X_train, y_train, rsID_train=None):
        """
        Train the model on the provided training data.

        Args:
            X_train: Training feature data.
            y_train: Training labels.
            rsID_train: Annotation list
        """
        self.model.fit(X_train, y_train)

    def evaluate(self, X_test, y_test, rsID_test=None):
        """
        Evaluate the model on test data and return performance metrics.

        Args:
            X_test: Test feature data.
            y_test: Test labels.
            rsID_test: Annotation list

        Returns:
            dict: Dictionary containing Mean Squared Error and R2 Score.
        """
        predictions = self.model.predict(X_test)
        # If predictions should be non-negative, set negatives to zero
        predictions = np.maximum(predictions, 0)
        mse = mean_squared_error(y_test, predictions)
        r2 = r2_score(y_test, predictions)
        return {'Mean Squared Error': mse, 'R2 Score': r2}

    def predict(self, X, y, rsID=None):

        predictions = self.model.predict(X)
        sorted_indices = np.argsort(predictions)
        ranks = np.argsort(sorted_indices) + 1
        if rsID is not None:
            rsID = np.array(rsID)
            # Create a dictionary mapping each rsID to its corresponding rank
            rsID_to_rank = dict(zip(rsID, ranks))
            # Sort rsIDs based on their ranks (ascending order)
            ordered_rsIDs = [x for _, x in sorted(zip(ranks, rsID), key=lambda pair: pair[0])]
            return predictions, ranks, rsID_to_rank, ordered_rsIDs
        return predictions, ranks
    def save_model(self, filepath):
        """
        Save the trained model to a file.

        Args:
            filepath (str): Path to the file where the model should be saved.
        """
        joblib.dump(self.model, filepath)
        print(f"Model saved to {filepath}")

    @classmethod
    def load_model(cls, filepath, config=None):
        """
        Load a model from a file.

        Args:
            filepath (str): Path to the file from which the model should be loaded.
            config (dict): Configuration dictionary for initializing the model structure.

        Returns:
            ClassicalModel: An instance of ClassicalModel with the loaded model.
        """
        # Initialize a new instance without setting self.model
        instance = cls(config or {})
        # Load the model state
        instance.model = joblib.load(filepath)
        print(f"Model loaded from {filepath}")
        return instance



# Model Factory to Load Models from Config
def load_model_from_config(config):
    method = config.get('method', 'neural_network')
    # seed = config.get('random_state')
    # seed_it(seed)
    if method == 'neural_network':
        return NeuralNetworkModel(config)
    else:
        return ClassicalModel(config)

# Example Main Function to Train and Evaluate the Model
def main(config):
    # Load Data
    X, y, rsID = load_data(r"C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\data")
    print(X.shape)


    # Split Data
    test_size = config.get('test_size', 0.2)
    random_state = config.get('random_state', 42)
    X_train, X_test, y_train, y_test, rsID_train, rsID_test = split_data(X, y, rsID, test_size)

    # Preprocess Data
    X_train_scaled, X_test_scaled, y_train_transformed, y_test_transformed, X_scaled, y_transformed = preprocess_data(X_train, X_test, y_train, y_test, X, y)
    # X_train_scaled, X_test_scaled = X_train, X_test
    print("X_train_scaled Mean:", np.mean(X_train_scaled), "Std:", np.std(X_train_scaled))
    print("X_test_scaled Mean:", np.mean(X_test_scaled), "Std:", np.std(X_test_scaled))
    print("y_train_transformed Mean:", np.mean(y_train_transformed), "Std:", np.std(y_train_transformed))
    print("y_test_transformed Mean:", np.mean(y_test_transformed), "Std:", np.std(y_test_transformed))
    # fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Plotting with adjusted axis limits to ignore extreme outliers
    # sns.histplot(X_train_scaled.flatten(), bins=300, kde=True, ax=axs[0, 0], color='skyblue', legend=False)
    # axs[0, 0].set_title('Distribution of X_train_scaled')
    # axs[0, 0].set_xlim([-5, 5])  # Adjusting x-axis limits
    #
    # sns.histplot(X_test_scaled.flatten(), bins=300, kde=True, ax=axs[0, 1], color='olive', legend=False)
    # axs[0, 1].set_title('Distribution of X_test_scaled')
    # axs[0, 1].set_xlim([-5, 5])  # Adjusting x-axis limits
    #
    # sns.histplot(y_train_transformed, bins=300, kde=True, ax=axs[1, 0], color='gold', legend=False)
    # axs[1, 0].set_title('Distribution of y_train_transformed')
    #
    # sns.histplot(y_test_transformed, bins=300, kde=True, ax=axs[1, 1], color='teal', legend=False)
    # axs[1, 1].set_title('Distribution of y_test_transformed')


    fig, axs = plt.subplots(1, 1, figsize=(12, 8))  # Only one row of plots

    # Histogram for X_scaled
    # sns.histplot(X_scaled.flatten(), bins=300, kde=True, color='skyblue', legend=False, ax=axs[0])
    # axs[0].set_title('Distribution of X_train_scaled')
    # axs[0].set_xlim([-5, 5])  # Adjusting x-axis limits

    # Histogram for y_transformed
    sns.histplot(y_transformed.flatten(), bins=300, kde=True, color='olive', legend=False, ax=axs[1])
    axs[1].set_title('Distribution of transformed score')

    plt.tight_layout()
    plt.show()

    plt.tight_layout()
    plt.show()
    # Load Model
    model = load_model_from_config(config)

    # Train Model
    method_name = config.get('method', 'neural_network').replace('_', ' ').title()
    print(f"\nTraining {method_name}...")
    model.train(X_train_scaled, y_train_transformed)

    # Evaluate Model
    print(f"\nEvaluating {method_name}...")
    results = model.evaluate(X_test_scaled, y_test_transformed)
    print(f"Evaluation Results: {results}")

    # retrieve the rank
    predictions, ranks, rsID_to_rank, ordered_rsIDs = model.predict(X_scaled, y_transformed, rsID)
    print(ordered_rsIDs[:200])
    names = rsid2gene(ordered_rsIDs[:200])
    print(names)


# 7. Example Configuration and Execution
if __name__ == "__main__":


    # Example config for Random Forest with hyperparameters
    config_rf = {
        'method': 'random_forest',
        'model_params': {
            'n_estimators': 100,
            'max_depth': None,
            'random_state': 42
        },
        'use_transformed_target_regressor': True,
        'test_size': 0.5,
        'random_state': 42
    }
    config_svr = {
        'method': 'svr',
        'test_size': 0.5,
        'random_state': 100
    }
    config_nn = {
        'method': 'neural_network',
        'test_size': 0.5,
        'random_state': 42
    }
    config_dt = {
        'method': 'decision_tree',
        'model_params': {
            'max_depth': 5,
            'random_state': 42
        },
        'use_transformed_target_regressor': True,
        'test_size': 0.5,
        'random_state': 42
    }

    config_lasso = {
        'method': 'lasso',
        'model_params': {
            'alpha': 0.1,
            'random_state': 42
        },
        'use_transformed_target_regressor': True,
        'test_size': 0.5,
        'random_state': 42
    }

    config_ridge = {
        'method': 'ridge',
        'model_params': {
            'alpha': 1.0,
            'random_state': 42
        },
        'use_transformed_target_regressor': True,
        'test_size': 0.5,
        'random_state': 42
    }

    config_lr= {
        'method': 'linear_regression',
        'model_params': {},
        'use_transformed_target_regressor': True,
        'test_size': 0.5,
        'random_state': 42
    }
    main(config_rf)
    # main(config_svr)
    # main(config_nn)
    # main(config_lr)
    # main(config_ridge)
    # main(config_lasso)
    # main(config_dt)
