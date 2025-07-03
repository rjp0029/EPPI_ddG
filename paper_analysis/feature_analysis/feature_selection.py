import pandas as pd
import numpy as np
from scipy.stats import entropy
# import a Linear Regressor, a Partial Least Squares Regressor, a Multi-layer Perceptron Regressor, a Support Vector Regressor, a K-Nearest-Neighbors Regressor, a Gradient Boosting Regressor, and a Random Forest Regressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.feature_selection import SelectFromModel, RFECV
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
import joblib
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

def recursive_feature_elimination(features, n_estimators=100, max_depth=30, min_samples_split=5, min_samples_leaf=2, max_features='sqrt'):
    # get the data in 20 fold splits making sure that the same mutations in different directions are in the same set
    # Create a function to define groups
    def extract_base_mutation(index):
        parts = index.split("_")
        pdb_id, mutation = parts[0], parts[1]
        reverse_mutation = f"{mutation[-1]}{mutation[1:-1]}{mutation[0]}"
        return f"{pdb_id}_{''.join(sorted([mutation, reverse_mutation]))}"
    
    X = features.drop(columns=['ddG'])
    y = features['ddG']

    groups = X.index.map(extract_base_mutation)
    group_splitter = GroupKFold(n_splits=20)
    model = RandomForestRegressor(random_state=RANDOM, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, max_features=max_features)
    rfecv = RFECV(estimator=model, step=1, cv=group_splitter.split(X, y, groups), scoring='r2', verbose=1, n_jobs=-1)

    rfecv.fit(X, y)
    print("Optimal number of features : %d" % rfecv.n_features_)

    # Plot number of features VS. cross-validation scores
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (r2)")
    plt.plot(range(1, len(rfecv.cv_results_['mean_test_score']) + 1), rfecv.cv_results_['mean_test_score'])
    

    # Get the selected features
    selected_features = X.columns[rfecv.support_]

    print("Selected features:", selected_features)

    # save the selected features to a file
    with open('selected_features.txt', 'w') as f:
        for feature in selected_features:
            f.write(f"{feature}\n")

    # Return the feature matrix with only the selected features
    return features[selected_features.tolist() + ['ddG']]  

if __name__ == '__main__':
    RANDOM = 167
    # file = '../data/single_tate_eppi_features.csv'
    file = '../data/ensemble_eppi_features.csv'
    eppi_features = load_features(file) 
    mutated_features = eppi_features
    print("shape of mutated features: ",len(mutated_features), ', ', len(mutated_features[list(mutated_features.keys())[0]]))
    file = '../data/wildtype_features.csv'
    eppi_features = load_features(file)
    original_features = eppi_features
    
    # get delta features (including the ddg from skempi)
    delta_features = get_delta_features(mutated_features, original_features)
    # print the shape of the delta features
    print('all delta features total', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # split the validation set off from the training set
    validation_set, delta_features = stratified_validation_split(delta_features)

    # permute mutations
    delta_features = permute_mutations(delta_features)
    print('all delta features permuted', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # reverse mutations
    delta_features = reverse_mutations(delta_features)
    print('all delta features permuted and reversed', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # reverse mutations
    validation_set = reverse_mutations(validation_set)
    print('validation set', len(validation_set), ', ', len(validation_set[list(validation_set.keys())[0]]))

    # convert to dataframe
    delta_features_df = pd.DataFrame(delta_features).T
    print(delta_features_df.shape)

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

    # Remove low entropy features
    delta_features_df = remove_low_entropy(delta_features_df, threshold=h)
    print(delta_features_df.shape, " after removing low entropy features")

    # Remove correlated features
    delta_features_df = remove_correlated_features(delta_features_df, threshold=r)
    print(delta_features_df.shape, " after removing correlated features")

    delta_features_df = recursive_feature_elimination(delta_features_df, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, max_features=max_features)
    print(delta_features_df.shape, " after recursive feature elimination")

    plt.show()

   


