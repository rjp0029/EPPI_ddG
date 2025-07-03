import joblib
import argparse
import warnings
warnings.filterwarnings("ignore")

def load_model(name="single_state_2421"):
    return joblib.load(f"PANTZ/models/{name}.pkl")

def load_delta_eppi_features():
    delta_eppi_features = []
    delta_feature_names = []
    with open("delta_eppi_features.txt", "r") as f:
        for line in f:
            delta_eppi_features.extend([float(line.split()[1])])
            delta_feature_names.append(line.split()[0][:-1])
    return delta_eppi_features, delta_feature_names

def check_feature_names(model_feature_names, delta_feature_names):
    if model_feature_names != delta_feature_names:
        raise ValueError("Feature names mismatch, terminating program.")

if __name__ == "__main__":
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="single_state_2421", help="Model name to load")
    args = parser.parse_args()
    # load model
    model = load_model(name=args.model)
    # load data
    delta_eppi_features, delta_feature_names = load_delta_eppi_features()
    # check if the feature names match
    check_feature_names(list(model.feature_names_in_), delta_feature_names)
    # run the prediction on the data
    prediction = model.predict([delta_eppi_features])
    print(prediction[0])
