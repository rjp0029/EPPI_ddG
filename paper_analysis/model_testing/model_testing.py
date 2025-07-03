import pandas as pd
import numpy as np
from scipy.stats import entropy
# import a Linear Regressor, a Partial Least Squares Regressor, a Multi-layer Perceptron Regressor, a Support Vector Regressor, a K-Nearest-Neighbors Regressor, a Gradient Boosting Regressor, and a Random Forest Regressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut, train_test_split
import re
sns.set_theme(style="whitegrid", context="talk", palette="muted")

# viridis color palette
viridis = sns.color_palette("viridis", as_cmap=True)
train_color = viridis(0.8)
test_color = viridis(0.33)

# validation_color = viridis(0.2)
validation_color = '#1F968B'
# flex_ddg_color = viridis(0.6)
flex_ddg_color = '#453781'

feature_names = {
    "n_epp_sb": "Number of persistent salt bridges",
    "n_epp_hb": "Number of persistent hydrogen bonds",
    "n_epp_hp": "Number of persistent hydrophobic interactions",
    "n_np_sb": "Number of non-persistent salt bridges",
    "n_np_hb": "Number of non-persistent hydrogen bonds",
    "n_np_hp": "Number of non-persistent hydrophobic interactions",
    "rl_epp_sb": "Rotamers lost in persistent salt bridges",
    "rl_epp_hb": "Rotamers lost in persistent hydrogen bonds",
    "rl_epp_hp": "Rotamers lost in persistent hydrophobic interactions",
    "rl_np_sb": "Rotamers lost in non-persistent salt bridges",
    "rl_np_hb": "Rotamers lost in non-persistent hydrogen bonds",
    "rl_np_hp": "Rotamers lost in non-persistent hydrophobic interactions",
    "dg_separated": "Rosetta Interaction Energy",
    "bsa": "Buried Surface Area",
    "sc": "Shape Complementarity",
    "epbe_sb": "Expected persistent binding energy of salt bridges",
    "epbe_hb": "Expected persistent binding energy of hydrogen bonds",
    "epbe_hp": "Expected persistent binding energy of hydrophobic interactions",
    "n_st_st_sb": "Number of stable-stable salt bridges",
    "n_st_st_hb": "Number of stable-stable hydrogen bonds",
    "n_st_st_hp": "Number of stable-stable hydrophobic interactions",
    "n_epp_polar": "Number of persistent polar interactions",
    "n_eppi": "Number of persistent interactions",
    "n_np_polar": "Number of non-persistent polar interactions",
    "n_nppi": "Number of non-persistent interactions",
    "rl_epp_polar": "Rotamers lost in persistent polar interactions",
    "rl_eppi": "Rotamers lost in persistent interactions",
    "rl_np_polar": "Rotamers lost in non-persistent polar interactions",
    "rl_nppi": "Rotamers lost in non-persistent interactions",
    "epbe_polar": "Expected persistent binding energy of polar interactions",
    "epbe": "Expected persistent binding energy",
    "npbe": "Non-persistent binding energy",
    "n_st_st_polar": "Number of stable-stable polar interactions",
    "n_st_st_pi": "Number of stable-stable interactions",
    "n_total": "Total number of interactions",
    "rl_total": "Total rotamers lost"
}

def load_skempi(file):
    # read in data from the skempi database
    df = pd.read_csv(file, sep=';')
    # drop rows where method isnt ITC or SPR or SFFL or SP or FL
    df = df[df['Method'].isin(['ITC', 'SPR', 'SFFL', 'SP', 'FL'])]
    # drop rows where mutation is not single point mutation
    df = df[df['Mutation(s)_PDB'].str.contains(',') == False]
    # convert kd (Molar) to kcal/mol using the equation: Kd = exp(dG/RT)
    R = 0.0019872041  # kcal/mol/K
    # convert temperature to int keeping only numeric values in the string
    df['Temperature'] = df['Temperature'].str.extract(r'(\d+)').astype(int)
    # convert kd to dG
    df['dG_wt'] = R * df['Temperature'] * np.log(df['Affinity_wt_parsed'])
    df['dG_mut'] = R * df['Temperature'] * np.log(df['Affinity_mut_parsed'])
    # calculate ddG of mutation
    df['ddG'] = df['dG_mut'] - df['dG_wt']
    # rank mutations by ddG
    df['rank'] = df['ddG'].rank()
    # sort by rank
    df = df.sort_values('rank')
    # remove rows with nan ddg
    df = df.dropna(subset=['ddG'])
    # drop those from bad pdbs
    # df = df[~df['#Pdb'].isin(FAILED_PREDICTION_PDBS + BAD_NUMBERING)]
    # set name to the first 4 characters of the pdb id and the mutation
    df['Name'] = df['#Pdb'] + "_" + df['Mutation(s)_PDB']
    # set the index to the pdb id and mutation
    df = df.set_index(['Name'])
    # remove the middle two parts of the Name, which are the interface chain names
    df.index = df.index.str.split('_').str[0] + "_" + df.index.str.split('_').str[-1]

    # # drop all mutations that have more than one entry
    # df = df[~df.index.isin(df.index.value_counts()[df.index.value_counts() > 1].index)]

    # order of methods for prioritization
    method_order = ['ITC', 'SPR', 'SFFL', 'SP', 'FL']
    
    # group by mutation and handle multiple entries
    def process_group(group):
        # sort the group by the method order
        group = group.sort_values(by='Method', key=lambda x: x.map({method: idx for idx, method in enumerate(method_order)}))
        # If there are multiple entries of the same method, take the average ddG
        group = group.groupby('Method', as_index=False).agg({'ddG': 'mean'}) 
        # sort again by the prioritized method order
        group = group.sort_values(by='Method', key=lambda x: x.map({method: idx for idx, method in enumerate(method_order)}))
        # Ttke the first method's ddG as the representative value
        return group.iloc[0]

    # Apply this processing to each mutation
    df = df.groupby(df.index).apply(process_group)
    # drop all columns except the index and ddG and method
    df = df[['ddG']]
    # sort the df on ddg values from highest to lowest
    df = df.sort_values('ddG', ascending=False)
    # print(df)
    return df

def load_features(file):
    # Load data as dictionary file looks like this:
    eppi_features = {}
    with open(file, 'r') as f:
        # the first line is the header which contains the keys
        keys = f.readline().strip().split(',')
        # remove the last key which is the end of the line
        keys = keys[:-1]
        for line in f:
            values = line.strip().split(',')
            # remove the last value which is the end of the line
            values = values[:-1]
            # create a dictionary with the keys and values
            eppi_features[values[0]] = {k: float(v) for k, v in zip(keys[1:], values[1:])}
    return eppi_features

def get_delta_features(mutated_features, original_features):
    # get the delta features
    delta_features = {}
    # go through the mutated features and get their change from the original features
    for entry in mutated_features:
        delta_features[entry] = {}
        for feature in mutated_features[entry]:
            delta_features[entry][feature] = mutated_features[entry][feature] - original_features[entry[:4]+"_filled"][feature]
        
    # load skempi
    skempi = load_skempi('../data/skempi_v2.csv')
    # get the ddg values for the mutations
    ddg = skempi.loc[delta_features.keys()]['ddG']
    # add the ddg values to the delta features
    for mutation in delta_features:
        delta_features[mutation]['ddG'] = ddg[mutation]

    return delta_features

def permute_mutations(delta_features):
    # Augment the delta features with inverse mutations and permutations
    augmented_delta_features = {}
    for mutation, delta_values_dict in delta_features.items():
        pdb = mutation.split('_')[0]
        position = mutation.split('_')[1][1:-1]
        # Calculate permutations with all other mutations
        for other_mutation, other_delta_values_dict in delta_features.items():
            other_pdb = other_mutation.split('_')[0]
            other_position = other_mutation.split('_')[1][1:-1]
            # skip if the mutation is not the same pdb and position
            if pdb != other_pdb or position != other_position:
                continue
            if mutation != other_mutation:
                # Direct permutation from mutation to other_mutation
                perm_key = other_pdb+"_"+mutation.split('_')[1][-1] + mutation.split('_')[1][1:-1] + other_mutation.split('_')[1][-1]
                augmented_delta_features[perm_key] = {}
                for feature, value in delta_values_dict.items():
                    augmented_delta_features[perm_key][feature] = other_delta_values_dict[feature] - value
    # add the augmented dictionary to the original dictionary
    delta_features.update(augmented_delta_features)
    return delta_features

def reverse_mutations(delta_features):
    # Augment the delta features with inverse mutations and permutations
    augmented_delta_features = {}
    print(len(delta_features))
    for mutation, delta_values_dict in delta_features.items():
        # Add inverse mutation
        reverse = mutation.split('_')[0] + '_' + mutation.split('_')[1][-1] + mutation.split('_')[1][1:-1] + mutation.split('_')[1][0]
        augmented_delta_features[reverse] = {}
        # add the reverse mutation to the dictionary
        for feature, value in delta_values_dict.items():
            augmented_delta_features[reverse][feature] = -value
    # add the augmented dictionary to the original dictionary
    delta_features.update(augmented_delta_features)
    return delta_features

def stratified_validation_split(delta_features, validation_size=0.2):
    validation_set = {}
    n_samples_validation = int(len(delta_features) * validation_size)

    # Get the number of mutations per pdb
    mutation_counts = {}
    for mutation in delta_features.keys():
        pdb_id = mutation.split('_')[0]
        if pdb_id not in mutation_counts:
            mutation_counts[pdb_id] = 0
        mutation_counts[pdb_id] += 1
    mutation_counts = {k: v for k, v in sorted(mutation_counts.items(), key=lambda item: item[1], reverse=True)}


    # go through each pdb with more than 10 mutations and sample 20%
    mutation_counts_greater = {k: v for k, v in mutation_counts.items() if v > 10}
    mutation_counts_less = {k: v for k, v in mutation_counts.items() if v <= 10}
    # print(f"Number of pdbs with more than 10 mutations: {len(mutation_counts_greater)}")
    # print(f"Number of pdbs with less than 10 mutations: {len(mutation_counts_less)}")

    np.random.seed(RANDOM)
    for pdb_id, count in mutation_counts_greater.items():
        n_samples_pdb = int(count * validation_size)
        pdb_mutations = [mutation for mutation in delta_features.keys() if mutation.split('_')[0] == pdb_id]
        samples = np.random.choice(pdb_mutations, n_samples_pdb, replace=False)
        for sample in samples:
            validation_set[sample] = delta_features[sample]

    # from the remaining pdbs sample 20% of the mutations
    n_samples_remaining = n_samples_validation - len(validation_set)
    remaining_mutations = [mutation for mutation in delta_features.keys() if mutation.split('_')[0] not in mutation_counts_greater]
    samples = np.random.choice(remaining_mutations, n_samples_remaining, replace=False)
    for sample in samples:
        validation_set[sample] = delta_features[sample]

    # # from the remaining pdbs sample 20% of the pdb complexes 
    # n_samples_remaining = n_samples_validation - len(validation_set)
    # remaining_pdbs = set([mutation.split('_')[0] for mutation in delta_features.keys()]) - set(mutation_counts_greater.keys())
    # samples = np.random.choice(list(remaining_pdbs), n_samples_remaining, replace=False)
    # for sample in samples:
    #     pdb_mutations = [mutation for mutation in delta_features.keys() if mutation.split('_')[0] == sample]
    #     sample = np.random.choice(pdb_mutations)
    #     validation_set[sample] = delta_features[sample]
    
    # Remove validation mutations from training
    for mutation in validation_set:
        del delta_features[mutation]

    training_mutations = delta_features.keys()
    training_positions = set([mutation.split('_')[0] + "_" + mutation.split('_')[1][1:-1] for mutation in delta_features.keys()])
    training_pdbs = set([mutation.split('_')[0] for mutation in delta_features.keys()])
    validation_mutations = validation_set.keys()
    validation_positions = set([mutation.split('_')[0] + "_" + mutation.split('_')[1][1:-1] for mutation in validation_set.keys()])
    validation_pdbs = set([mutation.split('_')[0] for mutation in validation_set.keys()])

    print("N unique mutations in the training set: ", len(training_mutations))
    print("N unique positions in the training set: ", len(training_positions))
    print("N unique pdbs in the training set: ", len(training_pdbs))

    print("N unique mutations in the validation set: ", len(validation_mutations))
    print("N unique positions in the validation set: ", len(validation_positions))
    print("N unique pdbs in the validation set: ", len(validation_pdbs))

    print("N unique mutations in training set that are not in validation set: ", len(set(training_mutations) - set(validation_mutations)))
    print("N unique positions in training set that are not in validation set: ", len(set(training_positions) - set(validation_positions)))
    print("N unique pdbs in training set that are not in validation set: ", len(set(training_pdbs) - set(validation_pdbs)))

    print("N unique mutations in validation set that are not in training set: ", len(set(validation_mutations) - set(training_mutations)))
    print("N unique positions in validation set that are not in training set: ", len(set(validation_positions) - set(training_positions)))
    print("N unique pdbs in validation set that are not in training set: ", len(set(validation_pdbs) - set(training_pdbs)))

    # plot the distributions of ddG values for the training and validation sets
    sns.kdeplot([delta_features[mutation]['ddG'] for mutation in delta_features], label='Training', color=train_color, fill=True, linewidth=2.5)
    sns.kdeplot([validation_set[mutation]['ddG'] for mutation in validation_set], label='Validation', color=validation_color, fill=True, linewidth=2.5)
    plt.xlabel('ΔΔG (kcal/mol)')
    plt.ylabel('Density')
    plt.legend()

    # # get the distribution of pdbs in the training and validation sets
    # training_pdb_counts = {mutation.split('_')[0]: 0 for mutation in delta_features.keys()}
    # for mutation in delta_features.keys():
    #     training_pdb_counts[mutation.split('_')[0]] += 1
    # # sort the dictionary by the number of mutations
    # training_pdb_counts = {k: v for k, v in sorted(training_pdb_counts.items(), key=lambda item: item[1], reverse=True)}

    # validation_pdb_counts = {mutation.split('_')[0]: 0 for mutation in validation_set.keys()}
    # for mutation in validation_set.keys():
    #     validation_pdb_counts[mutation.split('_')[0]] += 1

    # #plot the distribution of pdbs in the training and validation sets
    # plt.figure(figsize=(12, 6))
    # plt.bar(training_pdb_counts.keys(), training_pdb_counts.values(), label='Training')
    # plt.bar(validation_pdb_counts.keys(), validation_pdb_counts.values(), label='Validation')
    # plt.xlabel('PDB ID')
    # plt.ylabel('Number of mutations')
    # # rotate the x labels
    # plt.xticks(rotation=90)
    # plt.legend()
    return validation_set, delta_features

def remove_low_entropy(features, threshold=1):
    # remove low entropy features (< threshold) from the dataframe
    for feature in features.columns:
        if entropy(features[feature].value_counts(normalize=True), base=2) < threshold:
            features.drop(columns=[feature], inplace=True)
    return features

def remove_correlated_features(features, threshold=0.95):
    # remove correlated features (> threshold) from the dataframe
    corr_matrix = features.corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool_))
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
    features.drop(columns=to_drop, inplace=True)
    return features

def keep_selected_features(features, features_file):
    # load the features from the features file
    with open(features_file, 'r') as f:
        select_features = [line.strip() for line in f]
    # keep ddG
    select_features.append('ddG')
    # keep only the selected features
    features = features[select_features]
    return features

def balance_sample(delta_features, threshold=50, random_state=42):
    pdb_id_counts = delta_features.index.str.split("_").str[0].value_counts()
    
    pdb_id_counts_less = pdb_id_counts[pdb_id_counts <= 1000]
    pdb_id_counts_more = pdb_id_counts[pdb_id_counts > 1000]
    
    delta_features_balanced = delta_features[delta_features.index.str.split("_").str[0].isin(pdb_id_counts_less.index)].copy()
    
    def get_mutation_pair(index):
        parts = index.split("_")
        pdb_id, mutation = parts[0], parts[1]
        reverse_mutation = f"{mutation[-1]}{mutation[1:-1]}{mutation[0]}"
        return f"{pdb_id}_{''.join(sorted([mutation, reverse_mutation]))}"
    
    for pdb_id in pdb_id_counts_more.index:
        pdb_df = delta_features[delta_features.index.str.split("_").str[0] == pdb_id].copy()
        pdb_df["mutation_pair"] = pdb_df.index.map(get_mutation_pair)
        
        # Mean of absolute ddG per mutation pair
        pair_ddg = pdb_df.groupby("mutation_pair")["ddG"].apply(lambda x: np.mean(np.abs(x)))
        unique_pairs = pair_ddg.index.values
        
        # Weighting: Gaussian centered at 0 absolute ddG
        sigma = pair_ddg.std() if pair_ddg.std() > 0 else 1  # Avoid division by zero
        weights = np.exp(-0.5 * (pair_ddg.values / sigma) ** 2)
        # invert the weights
        weights = 1 - weights
        weights /= weights.sum()  # Normalize weights
        
        # Seed
        np.random.seed(random_state)
        sampled_pairs = np.random.choice(
            unique_pairs,
            size=min(threshold, len(unique_pairs)),
            replace=False,
            p=weights
        )
        
        pdb_df_sampled = pdb_df[pdb_df["mutation_pair"].isin(sampled_pairs)].drop(columns=["mutation_pair"])
        delta_features_balanced = pd.concat([delta_features_balanced, pdb_df_sampled])
    
    return delta_features_balanced

def direction_aware_kfold_split(X, y, n_splits = 10):
    # get the data in 20 fold splits making sure that the same mutations in different directions are in the same set
    # Create a function to define groups
    def extract_base_mutation(index):
        parts = index.split("_")
        pdb_id, mutation = parts[0], parts[1]
        reverse_mutation = f"{mutation[-1]}{mutation[1:-1]}{mutation[0]}"
        return f"{pdb_id}_{''.join(sorted([mutation, reverse_mutation]))}"
    
    # Extract group identifiers
    groups = X.index.map(extract_base_mutation)
    # Split the data
    group_splitter = GroupKFold(n_splits=n_splits)
    X_train_folds = []
    X_test_folds = []
    y_train_folds = []
    y_test_folds = []
    for train_index, test_index in group_splitter.split(X, y, groups=groups):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        X_train_folds.append(X_train)
        X_test_folds.append(X_test)
        y_train_folds.append(y_train)
        y_test_folds.append(y_test)
    return X_train_folds, X_test_folds, y_train_folds, y_test_folds

def position_aware_kfold_split(X, y, n_splits = 10):
    # get the data in group splits where one specific mutation position is left out each time
    # Create a function to define groups
    def extract_base_mutation(index):
        parts = index.split("_")
        pdb_id, mutation = parts[0], parts[1]
        return f"{pdb_id}_{mutation[1:-1]}"
    
    # Extract group identifiers
    groups = X.index.map(extract_base_mutation)
    # Split the data
    group_splitter = GroupKFold(n_splits=n_splits)
    X_train_folds = []
    X_test_folds = []
    y_train_folds = []
    y_test_folds = []
    for train_index, test_index in group_splitter.split(X, y, groups=groups):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        X_train_folds.append(X_train)
        X_test_folds.append(X_test)
        y_train_folds.append(y_train)
        y_test_folds.append(y_test)
    return X_train_folds, X_test_folds, y_train_folds, y_test_folds

def leave_one_complex_out_split(X, y):
    # get the data in group splits where one specific complex is left out each time
    # Create a function to define groups
    def extract_base_mutation(index):
        parts = index.split("_")
        pdb_id, mutation = parts[0], parts[1]
        return f"{pdb_id}"
    
    # Extract group identifiers
    groups = X.index.map(extract_base_mutation)
    # Split the data
    group_splitter = LeaveOneGroupOut()
    X_train_folds = []
    X_test_folds = []
    y_train_folds = []
    y_test_folds = []
    for train_index, test_index in group_splitter.split(X, y, groups=groups):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        X_train_folds.append(X_train)
        X_test_folds.append(X_test)
        y_train_folds.append(y_train)
        y_test_folds.append(y_test)
    return X_train_folds, X_test_folds, y_train_folds, y_test_folds

def train_and_test_model(X, y, n_estimators=100, max_depth=30, min_samples_split=5, min_samples_leaf=2, max_features='sqrt'):
    global k
    k=20
    # X_train_fold
    # X_train_folds, X_test_folds, y_train_folds, y_test_folds = direction_aware_kfold_split(X, y, n_splits=k)
    X_train_folds, X_test_folds, y_train_folds, y_test_folds = position_aware_kfold_split(X, y, n_splits=k)
    # X_train_folds, X_test_folds, y_train_folds, y_test_folds = leave_one_complex_out_split(X, y)

    # the best performing model has been the random forest regressor
    model = RandomForestRegressor(random_state=RANDOM, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, max_features=max_features)
    # model = RandomForestRegressor()
    # model = SVR()
    # model = LinearRegression()
    # model = GradientBoostingRegressor()
    # model = MLPRegressor()
    # model = KNeighborsRegressor()
    # model = PLSRegression()

    # # scale the data if using nn, svr, or knn
    # scaler = StandardScaler()
    # for i in range(len(X_train_folds)):
    #     X_train_folds[i] = scaler.fit_transform(X_train_folds[i])
    #     X_test_folds[i] = scaler.transform(X_test_folds[i])

    # store predictions from each fold
    kfold_y_pred = []
    kfold_y_pred_train = []
    kfold_y_test = []
    kfold_y_train = []

    feature_importances = {}

    for i in range(len(X_train_folds)):
        # get the training and testing data
        X_train = X_train_folds[i]
        X_test = X_test_folds[i]
        y_train = y_train_folds[i]
        y_test = y_test_folds[i]

        # print the sorted indexes of the training and testing data
        print(i+1, "train size: ", y_train.shape, y_train.shape[0] / (y_train.shape[0] + y_test.shape[0]), "test size: ", y_test.shape, y_test.shape[0] / (y_train.shape[0] + y_test.shape[0]))

        model.fit(X_train, y_train)
        
        # predict on train and test
        y_pred_train = model.predict(X_train)
        y_pred = model.predict(X_test)
        
        # extend the predictions
        kfold_y_pred.extend(y_pred)
        kfold_y_pred_train.extend(y_pred_train)
        kfold_y_test.extend(y_test)
        kfold_y_train.extend(y_train)

        # print top 10 feature importances
        feature_importances[i] = model.feature_importances_
        feature_importances[i] = {feature: importance for feature, importance in zip(X.columns, feature_importances[i])}
        feature_importances[i] = {k: v for k, v in sorted(feature_importances[i].items(), key=lambda item: item[1], reverse=True)}

    
    # plot the predicted vs. actual ddG
    y_pred = np.array(kfold_y_pred)
    y_pred_train = np.array(kfold_y_pred_train)
    y_test = np.array(kfold_y_test)
    y_train = np.array(kfold_y_train)

    return y_pred, y_pred_train, y_test, y_train, model, feature_importances

def get_fc(prediction, experimental):
    # fraction correct gets the number fraction of mutations that are correctly desgtabilizing (less than -1) neutral (-1 to 1) and stabilizing (greater than 1)
    # get the number of mutations that are correctly destabilizing, neutral and stabilizing
    correct_destabilizing = np.sum((experimental < -1) & (prediction < -1))
    correct_neutral = np.sum((experimental >= -1) & (experimental <= 1) & (prediction >= -1) & (prediction <= 1))
    correct_stabilizing = np.sum((experimental > 1) & (prediction > 1))
    # get the total number of mutations
    total = len(experimental)
    # get the fraction correct
    fc = (correct_destabilizing + correct_neutral + correct_stabilizing) / total
    return fc

def plot_model(y_pred, y_pred_train, y_test, y_train, model, feature_importances):
    # calculate training metrics
    train_mse = mean_squared_error(y_train, y_pred_train)
    train_r_squared = 1 - (np.sum((y_train - y_pred_train) ** 2) / np.sum((y_train - np.mean(y_train)) ** 2))
    train_pearson = np.corrcoef(y_train, y_pred_train)[0, 1]
    train_mae = np.mean(np.abs(y_train - y_pred_train))

    # get the test metrics
    test_mse = mean_squared_error(y_test, y_pred)
    test_r_squared = 1 - (np.sum((y_test - y_pred) ** 2) / np.sum((y_test - np.mean(y_test)) ** 2))
    test_pearson = np.corrcoef(y_test, y_pred)[0, 1]
    test_mae = np.mean(np.abs(y_test - y_pred))

    # print metrics
    print("Training:", "MSE:", train_mse, "R^2:", train_r_squared, "Pearson:", train_pearson, "MAE:", train_mae)
    print(model.__class__.__name__, "Testing:", "MSE:", test_mse, "R^2:", test_r_squared, "Pearson:", test_pearson, "MAE:", test_mae)

    # print the predictions correct within their respective thresholds
    within_point4_kcal = np.sum(np.abs(y_test - y_pred) <= 0.4) / len(y_test)
    within_1_kcal = np.sum(np.abs(y_test - y_pred) <= 1) / len(y_test)
    correct_sign = np.sum(np.sign(y_test) == np.sign(y_pred)) / len(y_test)
    print(f"Within 0.4 kcal/mol: {within_point4_kcal:.2%}")
    print(f"Within 1 kcal/mol: {within_1_kcal:.2%}")
    print(f"Correct sign: {correct_sign:.2%}")

    # get bands for different KD fold changes
    R = 1.9872036e-3
    RT = R * 298.15
    kd_fold_changes = [1, 5, 10, 50, 100, 500, 1000]
    thresholds = [RT * np.log(fc) for fc in kd_fold_changes]
    sign_accuracies = []
    counts = []

    print("\nSign accuracy for KD thresholds:")
    for fc, threshold in zip(kd_fold_changes, thresholds):
        # print(fc, threshold)
        mask = np.abs(y_test) >= threshold
        count = np.sum(mask)
        if count == 0:
            print(f"No mutations with experimental ΔΔG > {threshold:.3f} kcal/mol ({fc}x Kd change)")
            sign_accuracies.append(np.nan)
            counts.append(0)
            continue
        # print(count)
        correct_sign_fraction = np.sum(np.sign(y_test[mask]) == np.sign(y_pred[mask])) / count
        # print(correct_sign_fraction)
        sign_accuracies.append(correct_sign_fraction)
        counts.append(count)
        print(f"  > {threshold:.3f} kcal/mol ({fc}x Kd change): {correct_sign_fraction:.2%} ({count} mutations)")

    # bar chart of sign accuracies for different thresholds
    plt.figure(figsize=(8, 5))
    sns.barplot(
        x=[f"{fc}x" for fc in kd_fold_changes],
        y=sign_accuracies,
        color="grey")
    plt.ylim(0, 1.1)
    plt.ylabel("Sign Accuracy")
    plt.xlabel("Kd Fold Change Threshold")
    plt.title("Sign Accuracy for Different Kd Fold Change Thresholds")
    # for i, (acc, count) in enumerate(zip(sign_accuracies, counts)):
    #     if not np.isnan(acc):
    #         plt.text(i, acc + 0.02, f"{acc:.1%}\n(n={count})", ha='center')

    min_val = min(y_test.min(), y_pred.min(), y_train.min(), y_pred_train.min())
    max_val = max(y_test.max(), y_pred.max(), y_train.max(), y_pred_train.max())
    plt.figure(figsize=(10, 8))
    # ideal line
    plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='black', label='y = x')
    # scatter plots
    sns.scatterplot(x=y_train, y=y_pred_train, color=train_color, s=50)
    sns.scatterplot(x=y_test, y=y_pred, color=test_color, s=60)
    # # plot the train fit line
    # m, b = np.polyfit(y_train, y_pred_train, 1)
    # plt.plot([min_val, max_val], [m * min_val + b, m * max_val + b], color=fit_line_color, label=f'Training Fit (y = {m:.2f}x + {b:.2f})')
    # put lines at x=0 and y=0
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)
    # set the tick labels to be 0, 2, 4, 6, 8, 10
    plt.xticks(np.arange(-10, 11, 5), labels=[f"{i}" for i in np.arange(-10, 11, 5)], fontsize=28)
    plt.yticks(np.arange(-10, 11, 5), labels=[f"{i}" for i in np.arange(-10, 11, 5)], fontsize=28)
    
    # set the axis labels
    plt.xlabel('Experimental ΔΔG (kcal/mol)', fontsize=24)
    plt.ylabel('Predicted ΔΔG (kcal/mol)', fontsize=24)
    # # set the title
    # plt.title(
    #     f'{model.__class__.__name__} ΔΔG Prediction\nTest Pearson R: {test_pearson:.3f} | Test $R^2$: {test_r_squared:.3f}',
    #     fontsize=19,
    #     pad=25)
    # set the legend
    plt.legend(['y = x', 
                f'Training (R={train_pearson:.2f})',
                f'Testing (R={test_pearson:.2f})'],
                frameon=True,
                fontsize=20,
                title_fontsize=20,
                loc='upper left')

    plt.grid(True, which='both', linestyle='--', linewidth=0.3)
    plt.tight_layout()
    plt.savefig('../figures/model_performance.png', dpi=1200, bbox_inches='tight')

    # # plot the top 15 feature importances
    # plt.figure(figsize=(10, 8))
    # sns.barplot(x=list(feature_importances[0].values())[:15], y=list(feature_importances[0].keys())[:15], color='grey')
    # plt.xlabel('Feature Importance')
    # plt.ylabel('Feature')
    # plt.title('Top 15 Feature Importances')
    # plt.tight_layout()

    # average the feature importances
    feature_importances = {feature: np.mean([importances[feature] for importances in feature_importances.values()]) for feature in feature_importances[0].keys()}
    # sort the feature importances
    feature_importances = {k: v for k, v in sorted(feature_importances.items(), key=lambda item: item[1], reverse=True)}

    # write the feature importances to a csv file
    with open(model.__class__.__name__ + '_feature_importances.txt', 'w') as f:
        for feature, importance in feature_importances.items():
            parts = re.split(r'([*/+-])', feature)
            readable_parts = [feature_names.get(part, part) for part in parts]
            readable_feature = "".join(readable_parts)
            f.write(f"{readable_feature}: {importance:.4f}\n")

    return test_pearson, test_mae, get_fc(y_pred, y_test)

def validate_model(model, X_train, y_train, X_validate, y_validate):

    # Verify that reversed mutations are properly paired
    for mutation in X_validate.index:
        pdb = mutation.split('_')[0]
        mut_str = mutation.split('_')[1]
        reverse = f"{pdb}_{mut_str[-1]}{mut_str[1:-1]}{mut_str[0]}"
        if reverse not in X_validate.index:
            print(f"Warning: Missing reverse mutation for {mutation}")

    # # Save the indices of the training and validation data
    # X_train_indices = X_train.index
    # X_validate_indices = X_validate.index

    # retrain the model on all the training data
    model.fit(X_train, y_train)
    # predict on the training data
    y_pred_train = model.predict(X_train)
    negative_X_validate = -X_validate
    # predict on the validation data
    y_pred_validate = model.predict(X_validate)
    # predict on the negative validation data
    y_pred_negative_validate = model.predict(negative_X_validate)
    # # average the predictions
    # y_pred_validate = (y_pred_validate - y_pred_negative_validate) / 2

    # get performace on pdb codes from the validation set that are seen in the training set
    # get the pdbs in the training set
    training_pdbs = set([mutation.split('_')[0] for mutation in X_train.index])
    # get the pdbs in the validation set
    validation_pdbs = set([mutation.split('_')[0] for mutation in X_validate.index])
    # get the pdbs in the validation set that are in the training set
    common_pdbs = training_pdbs.intersection(validation_pdbs)
    # get the indices of the common pdbs
    indices = [i for i, mutation in enumerate(X_validate.index) if mutation.split('_')[0] in common_pdbs]
    # get the predictions for the common pdbs
    y_pred_validate_common = y_pred_validate[indices]
    y_validate_common = y_validate.iloc[indices]
    # Calculate metrics
    test_r_squared = 1 - (np.sum((y_validate_common - y_pred_validate_common) ** 2) / np.sum((y_validate_common - np.mean(y_validate_common)) ** 2))
    test_pearson = np.corrcoef(y_validate_common, y_pred_validate_common)[0, 1]
    test_mae = np.mean(np.abs(y_validate_common - y_pred_validate_common))
    test_fc = get_fc(y_pred_validate_common, y_validate_common)
    print("Number of mutations from common pdbs: ", len(y_validate_common))
    print("Seen PDBs Pearson: ", test_pearson, "R^2: ", test_r_squared, "MAE: ", test_mae, "FC: ", test_fc)

    # get performace on pdb codes from the validation set that are not seen in the training set
    # get the pdbs in the validation set that are not in the training set
    unseen_pdbs = validation_pdbs - training_pdbs
    # get the indices of the unseen pdbs
    indices = [i for i, mutation in enumerate(X_validate.index) if mutation.split('_')[0] in unseen_pdbs]
    # get the predictions for the unseen pdbs
    y_pred_validate_unseen = y_pred_validate[indices]
    y_validate_unseen = y_validate.iloc[indices]
    # Calculate metrics
    test_r_squared = 1 - (np.sum((y_validate_unseen - y_pred_validate_unseen) ** 2) / np.sum((y_validate_unseen - np.mean(y_validate_unseen)) ** 2))
    test_pearson = np.corrcoef(y_validate_unseen, y_pred_validate_unseen)[0, 1]
    test_mae = np.mean(np.abs(y_validate_unseen - y_pred_validate_unseen))
    test_fc = get_fc(y_pred_validate_unseen, y_validate_unseen)
    print("Number of mutations from unseen pdbs: ", len(y_validate_unseen))
    print("\n\nUnseen PDBs Pearson: ", test_pearson, "R^2: ", test_r_squared, "MAE: ", test_mae, "FC: ", test_fc, "\n\n")    
    return model, y_pred_train, y_pred_validate, test_pearson, test_r_squared

def validate_model_positions(model, X_train, y_train, X_validate, y_validate):
    # get set of positions (PDB_position) from training data, ignoring mutation direction
    training_positions = set([f"{mutation.split('_')[0]}_{mutation.split('_')[1][2:-1]}" 
                            for mutation in X_train.index])
    
    # get set of positions from validation data, ignoring mutation direction
    validation_positions = set([f"{mutation.split('_')[0]}_{mutation.split('_')[1][2:-1]}"
                              for mutation in X_validate.index])

    # get positions in validation set that were seen in training
    common_positions = training_positions.intersection(validation_positions)
    unseen_positions = validation_positions - training_positions
    
    # # print sorted positions with counts
    # print("\nTraining positions:")
    # training_pos_counts = {}
    # for mutation in X_train.index:
    #     pos = f"{mutation.split('_')[0]}_{mutation.split('_')[1][2:-1]}"
    #     training_pos_counts[pos] = training_pos_counts.get(pos, 0) + 1
    # for pos in sorted(training_positions):
    #     print(f"{pos}: {training_pos_counts[pos]} mutations")
        
    # print("\nValidation positions seen in training:")
    # validation_pos_counts = {}
    # for mutation in X_validate.index:
    #     pos = f"{mutation.split('_')[0]}_{mutation.split('_')[1][2:-1]}"
    #     validation_pos_counts[pos] = validation_pos_counts.get(pos, 0) + 1
    # for pos in sorted(common_positions):
    #     print(f"{pos}: {validation_pos_counts[pos]} mutations")
        
    # print("\nValidation positions not seen in training:")
    # for pos in sorted(unseen_positions):
    #     print(f"{pos}: {validation_pos_counts[pos]} mutations")

    # get predictions for validation set
    y_pred_validate = model.predict(X_validate)
    
    # get indices for mutations at common positions
    common_indices = [i for i, mutation in enumerate(X_validate.index) 
                     if f"{mutation.split('_')[0]}_{mutation.split('_')[1][2:-1]}" in common_positions]
    
    # calculate metrics for seen positions
    y_pred_validate_common = y_pred_validate[common_indices]
    y_validate_common = y_validate.iloc[common_indices]
    
    seen_r_squared = 1 - (np.sum((y_validate_common - y_pred_validate_common) ** 2) / 
                         np.sum((y_validate_common - np.mean(y_validate_common)) ** 2))
    seen_pearson = np.corrcoef(y_validate_common, y_pred_validate_common)[0, 1]
    seen_mae = np.mean(np.abs(y_validate_common - y_pred_validate_common))
    seen_fc = get_fc(y_pred_validate_common, y_validate_common)
    
    print("\nSeen Positions Performance:")
    print(f"Number of mutations from seen positions: {len(y_validate_common)}")
    print(f"Pearson: {seen_pearson:.5f}, R^2: {seen_r_squared:.5f}, MAE: {seen_mae:.5f}, FC: {seen_fc:.5f}")

    # get indices for mutations at unseen positions
    unseen_indices = [i for i, mutation in enumerate(X_validate.index)
                     if f"{mutation.split('_')[0]}_{mutation.split('_')[1][2:-1]}" in unseen_positions]
    
    # calculate metrics for unseen positions
    y_pred_validate_unseen = y_pred_validate[unseen_indices]
    y_validate_unseen = y_validate.iloc[unseen_indices]
    
    unseen_r_squared = 1 - (np.sum((y_validate_unseen - y_pred_validate_unseen) ** 2) /
                           np.sum((y_validate_unseen - np.mean(y_validate_unseen)) ** 2))
    unseen_pearson = np.corrcoef(y_validate_unseen, y_pred_validate_unseen)[0, 1]
    unseen_mae = np.mean(np.abs(y_validate_unseen - y_pred_validate_unseen))
    unseen_fc = get_fc(y_pred_validate_unseen, y_validate_unseen)
    
    print("\nUnseen Positions Performance:") 
    print(f"Number of mutations from unseen positions: {len(y_validate_unseen)}")
    print(f"Pearson: {unseen_pearson:.5f}, R^2: {unseen_r_squared:.5f}, MAE: {unseen_mae:.5f}, FC: {unseen_fc:.5f}")

    return (seen_pearson, seen_mae, seen_fc), (unseen_pearson, unseen_mae, unseen_fc)

def compare_flex_ddg_validation(y_pred_validate, y_validate):
    # read in the flex_ddgs:
    flex_ddgs = pd.read_csv('../data/flex_ddgs.csv', index_col=1)

    # print the total correlation for the flex_ddgs
    print("Flex ddG total Pearson: ", np.corrcoef(flex_ddgs['total_score'], flex_ddgs['ddG'])[0, 1])

    # get the reverse mutations for flex_ddgs
    for mut in flex_ddgs.index:
        pdb_id = mut.split('_')[0]
        mutation = mut.split('_')[1]
        reverse_mutation = f"{mutation[-1]}{mutation[1:-1]}{mutation[0]}"
        # add the reverse mutation to the flex_ddgs and its negative ddg and total score
        flex_ddgs.loc[f"{pdb_id}_{reverse_mutation}"] = [0, -flex_ddgs.loc[mut]['total_score'], -flex_ddgs.loc[mut]['ddG']]

    # print the total correlation for the flex_ddgs
    print("Reversed Flex ddG total Pearson: ", np.corrcoef(flex_ddgs['total_score'], flex_ddgs['ddG'])[0, 1])    

    # get the validation and flex_ddg in the same order
    flex_ddgs = flex_ddgs.loc[y_validate.index]

    # get the mutations in both the flex_ddgs and the validation set
    flex_ddg_corr = np.corrcoef(flex_ddgs['total_score'], flex_ddgs['ddG'])[0, 1]
    flex_ddg_r_squared = 1 - (np.sum((flex_ddgs['total_score'] - flex_ddgs['ddG']) ** 2) / np.sum((flex_ddgs['total_score'] - np.mean(flex_ddgs['total_score'])) ** 2))
    flex_ddg_mae = np.mean(np.abs(flex_ddgs['total_score'] - flex_ddgs['ddG']))
    flex_ddg_fc = get_fc(flex_ddgs['total_score'], flex_ddgs['ddG'])
    print("Flex ddG Pearson: ", flex_ddg_corr, "R^2: ", flex_ddg_r_squared, "MAE: ", flex_ddg_mae, "FC: ", flex_ddg_fc)

    # Calculate metrics
    test_r_squared = 1 - (np.sum((y_validate - y_pred_validate) ** 2) / np.sum((y_validate - np.mean(y_validate)) ** 2))
    test_pearson = np.corrcoef(y_validate, y_pred_validate)[0, 1]
    test_mae = np.mean(np.abs(y_validate - y_pred_validate))
    test_fc = get_fc(y_pred_validate, y_validate)
    print("Validation Pearson: ", test_pearson, "R^2: ", test_r_squared, "MAE: ", test_mae, "FC: ", test_fc)

    #plot the flex_ddg vs the validation data 
    plt.figure(figsize=(10, 8))
    min_val = min(y_validate.min(), y_pred_validate.min(), flex_ddgs['ddG'].min(), flex_ddgs['total_score'].min())
    max_val = max(y_validate.max(), y_pred_validate.max(), flex_ddgs['ddG'].max(), flex_ddgs['total_score'].max())
    # ideal line
    plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='black', label='y = x')
    sns.scatterplot(x=flex_ddgs['ddG'], y=flex_ddgs['total_score'], color=flex_ddg_color, s=50)
    sns.scatterplot(x=y_validate, y=y_pred_validate, color=validation_color, s=60)
    # set the tick labels to be 0, 2, 4, 6, 8, 10
    plt.xticks(np.arange(-10, 11, 5), labels=[f"{i}" for i in np.arange(-10, 11, 5)], fontsize=28)
    plt.yticks(np.arange(-10, 11, 5), labels=[f"{i}" for i in np.arange(-10, 11, 5)], fontsize=28)
    # put lines at x=0 and y=0
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)
    plt.legend(['y = x',
                f'flex_ddg (R={flex_ddg_corr:.2f})',
                f'EPPI_ddG (R={test_pearson:.2f})'],
                frameon=True,
                fontsize=20,
                title_fontsize=20,
                loc='upper left')
    plt.xlabel('Experimental ΔΔG (kcal/mol)', fontsize=24)
    plt.ylabel('Predicted ΔΔG (kcal/mol)', fontsize=24)
    plt.grid(True, which='both', linestyle='--', linewidth=0.3)
    plt.tight_layout()
    # plt.title('EPPI_ddg vs. flex_ddg Validation ΔΔG Prediction', fontsize=19, pad=25)
    plt.savefig('../figures/flex_ddg_vs_eppi_ddg.png', dpi=1200, bbox_inches='tight')

    return test_pearson - flex_ddg_corr

if __name__ == '__main__':
    RANDOM = 167
    # file = '../data/single_state_eppi_features.csv'
    file = '../data/ensemble_eppi_features.csv'
    eppi_features = load_features(file)
    mutated_features = eppi_features
    print("shape of mutated features: ",len(mutated_features), ', ', len(mutated_features[list(mutated_features.keys())[0]]))
    file = '../data/wildtype_features.csv'
    eppi_features = load_features(file)
    original_features = eppi_features
    
    # get delta features (including the ddg from skempi)
    delta_features = get_delta_features(mutated_features, original_features)
    # print(delta_features.keys())
    # delta_features = {k: v for k, v in delta_features.items() if '1CHO' not in k}    # print the shape of the delta features
    print('all delta features total', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # split the validation set off from the training set
    validation_set, delta_features = stratified_validation_split(delta_features)

    # delta_features_list = list(delta_features.keys())
    # delta_features_list = [x[:4] for x in delta_features_list]
    # delta_features_pdbs = list(set(delta_features_list))
    # print(delta_features_pdbs)

    # validation_list = list(validation_set.keys())
    # validation_list = [x[:4] for x in validation_list]
    # validation_pdbs = list(set(validation_list))
    # print(validation_pdbs)

    # # Convert to sets for easy set operations
    # train_set = set(delta_features_pdbs)
    # val_set = set(validation_pdbs)

    # # Compute intersections and differences
    # overlap = train_set & val_set
    # train_only = train_set - val_set
    # val_only = val_set - train_set

    # # Print the counts
    # print(f'Number of PDBs only in training set: {len(train_only)}')
    # print(f'Number of PDBs only in validation set: {len(val_only)}')
    # print(f'Number of PDBs in both training and validation sets: {len(overlap)}')


    # # write the validation mutations (just indexes) to a file
    # validation_mutations = list(validation_set.keys())
    # # sort the mutations
    # validation_mutations.sort()
    # with open('validation_mutations_set.txt', 'w') as f:
    #     for i, mutation in enumerate(validation_mutations):
    #         # write the mutation and the ddg value
    #         f.write(f"{mutation}, {validation_set[mutation]['ddG']:.2f}\n")


    # # write the training mutations (just indexes) to a file
    # delta_mutations = list(delta_features.keys())
    # # sort the mutations
    # delta_mutations.sort()
    # with open('delta_features_set.txt', 'w') as f:
    #     for i, mutation in enumerate(delta_mutations):
    #         # write the mutation and the ddg value
    #         f.write(f"{mutation}, {delta_features[mutation]['ddG']:.2f}\n")


    # permute mutations
    delta_features = permute_mutations(delta_features)
    print('all delta features permuted', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # reverse mutations
    delta_features = reverse_mutations(delta_features)
    print('all delta features permuted and reversed', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # # reverse mutations
    # validation_set = reverse_mutations(validation_set)
    # print('validation set', len(validation_set), ', ', len(validation_set[list(validation_set.keys())[0]]))

    # convert to dataframe
    delta_features_df = pd.DataFrame(delta_features).T
    print(delta_features_df.shape)
    # # write to csv  
    # delta_features_df.to_csv('all_delta_features.csv')

    # # load the data
    # delta_features_df = pd.read_csv('all_delta_features.csv', index_col=0)

    # HYPERPARAMETERS
    # parameters for the data processing
    h = 1
    r = 0.90
    b = 50
    # hyperparameters for the random forest regressor
    n_estimators = 300
    max_depth = 30
    min_samples_split = 2
    min_samples_leaf = 1
    max_features = 'sqrt'

    # # Remove low entropy features
    # delta_features_df = remove_low_entropy(delta_features_df, threshold=h)
    # print(delta_features_df.shape, " after removing low entropy features")

    # # Remove correlated features
    # delta_features_df = remove_correlated_features(delta_features_df, threshold=r)
    # print(delta_features_df.shape, " after removing correlated features")

    # keep only select features from recursive feature elimination
    delta_features_df = keep_selected_features(delta_features_df, '../data/selected_features.txt')
    print(delta_features_df.shape, " after keeping select features")

    # resample overrepresented complexes to balance the dataset
    delta_features_df = balance_sample(delta_features_df, threshold=b)
    print(delta_features_df.shape, " after balancing the data")

    # # write the balanced data to a file
    # delta_features_df.to_csv('selected_delta_features.csv')
    # delta_features_df = pd.read_csv('selected_delta_features.csv', index_col=0)

    # get X and y
    y = delta_features_df['ddG']
    X = delta_features_df.drop(columns=['ddG'])

    # train and test model
    y_pred, y_pred_train, y_test, y_train, model, feature_importances = train_and_test_model(X, y, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, max_features=max_features)

    # Evaluate model
    test_pearson_train, test_mae_train, fc_train = plot_model(y_pred, y_pred_train, y_test, y_train, model, feature_importances)

    # # get the validation datas df
    # validation_delta_features_df = pd.DataFrame(validation_set).T
    # # keep only the same fetures that are in the training data (and ddG)
    # validation_delta_features_df = validation_delta_features_df[X.columns.tolist() + ['ddG']]
    # # get the validation data
    # y_validate = validation_delta_features_df['ddG']
    # X_validate = validation_delta_features_df.drop(columns=['ddG'])

    # # validate the model
    # model, y_pred_train, y_pred_validate, unseen_pearson, unseen_r_squared = validate_model(model, X, y, X_validate, y_validate)

    # # also validate based on positions
    # seen_metrics, unseen_metrics = validate_model_positions(model, X, y, X_validate, y_validate)

    # print(seen_metrics)
    # print(unseen_metrics)

    # # evaluate the validation data
    # test_pearson_validation, test_mae_validation, fc_validation = plot_model(y_pred_validate, y_pred_train, y_validate, y, model, feature_importances)
    # delta_r = compare_flex_ddg_validation(y_pred_validate, y_validate)
    # print("Validation Delta R: ", delta_r)
    # # plt.show()
