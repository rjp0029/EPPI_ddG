import pandas as pd
import numpy as np
from scipy.stats import entropy
# import a Linear Regressor, a Partial Least Squares Regressor, a Multi-layer Perceptron Regressor, a Support Vector Regressor, a K-Nearest-Neighbors Regressor, a Gradient Boosting Regressor, and a Random Forest Regressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.decomposition import PCA
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut, train_test_split
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import euclidean_distances
import umap
from sklearn.manifold import TSNE
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
    print(pdb_id_counts)
    
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

    # get the unselected mutations
    unselected_mutations = delta_features[~delta_features.index.isin(delta_features_balanced.index)]
    print(f"Number of unselected mutations: {len(unselected_mutations)}")
    
    return delta_features_balanced

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

def test_leakage(X, y, n_estimators=100, max_depth=30, min_samples_split=5, min_samples_leaf=2, max_features='sqrt'):
    global k
    k = 20

    X_train_folds, X_test_folds, y_train_folds, y_test_folds = leave_one_complex_out_split(X, y)

    model = RandomForestRegressor(
        random_state=RANDOM, 
        n_estimators=n_estimators, 
        max_depth=max_depth, 
        min_samples_split=min_samples_split, 
        min_samples_leaf=min_samples_leaf, 
        max_features=max_features
    )

    # Store predictions & performance metrics
    kfold_y_pred, kfold_y_pred_train = [], []
    kfold_y_test, kfold_y_train = [], []
    pccs, pcc_improvements = {}, {}
    feature_importances = {}
    leakage_tracking = {}

    # Sort folds by training set size for debugging purposes
    sorted_folds = sorted(zip(X_train_folds, X_test_folds, y_train_folds, y_test_folds), key=lambda x: x[0].shape[0])
    X_train_folds, X_test_folds, y_train_folds, y_test_folds = zip(*sorted_folds)

    print("Number of test samples in each fold:", [X_test.shape[0] for X_test in X_test_folds])

    improved_pdbs = []
    worse_pdbs = []
    n_leaked_pairs = 6

    for i in range(len(X_train_folds)):
        X_train, X_test = X_train_folds[i], X_test_folds[i]
        y_train, y_test = y_train_folds[i], y_test_folds[i]

        print(f"{i+1} Pre-leak train size: {y_train.shape}, {y_train.shape[0] / (y_train.shape[0] + y_test.shape[0])}")
        print(f"Test size: {y_test.shape}, {y_test.shape[0] / (y_train.shape[0] + y_test.shape[0])}")

        # Scale the features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_test)
        X_scaled = pd.DataFrame(X_scaled, index=X_test.index, columns=X_test.columns)

        # Perform KMeans clustering on the test set to select the best representative mutations
        num_clusters = min(n_leaked_pairs, len(X_scaled))  # Select 5 clusters or fewer if less data is available
        kmeans = KMeans(n_clusters=num_clusters, random_state=RANDOM).fit(X_scaled)

        # Find the closest mutation to each cluster center
        cluster_centers = kmeans.cluster_centers_
        cluster_representatives = [
            X_test.iloc[np.argmin(np.linalg.norm(X_scaled - center, axis=1))].name 
            for center in cluster_centers
        ]

        # Select 10 representative mutation pairs
        representative_mutations = cluster_representatives[:n_leaked_pairs]

        # Identify reverse mutations and remove both forward and reverse mutations from the test set
        mutations_to_remove = set(representative_mutations)
        for mutation in representative_mutations:
            parts = mutation.split("_")
            complex_id, mutation_info = parts[0], parts[1]
            rev_mutation = mutation_info[0] + mutation_info[-1] + mutation_info[1:-1]
            reverse_mutation_id = f"{complex_id}_{rev_mutation}"
            if reverse_mutation_id in X_test.index:
                mutations_to_remove.add(reverse_mutation_id)

        # Create validation set by removing representative mutations from test set
        X_validation = X_test.drop(list(mutations_to_remove))
        y_validation = y_test.drop(list(mutations_to_remove))

        # Test the model on the validation set before leaking any pairs
        model.fit(X_train, y_train)
        y_pred = model.predict(X_validation)

        # Compute PCC for pre-leak validation set
        pcc_pre = np.corrcoef(y_validation, y_pred)[0, 1]
        pdb_id = X_validation.index[0].split("_")[0]
        pccs[pdb_id + "_preleak"] = pcc_pre
        print(f"{pdb_id} Pre-leak PCC: {pcc_pre}")

        # Track PCC and MAE for incremental leakage
        leakage_tracking[pdb_id] = {'pcc': [pcc_pre], 'mae': [np.mean(np.abs(y_pred - y_validation))], 'leaked_pairs': [0]}

        pairs_leaked = 0

        while pairs_leaked < n_leaked_pairs:
            # Select the next mutation to leak from the representative mutations
            next_mutation = representative_mutations[pairs_leaked]

            # Ensure reverse mutation is included
            parts = next_mutation.split("_")
            complex_id, mutation_info = parts[0], parts[1]
            
            rev_mutation = mutation_info[0] + mutation_info[-1] + mutation_info[1:-1]  
            reverse_mutation_id = f"{complex_id}_{rev_mutation}"

            # Add both mutations if the reverse exists in test set
            mutation_pair = {next_mutation, reverse_mutation_id} if reverse_mutation_id in X_test.index else {next_mutation}

            # Add selected mutations to training set
            X_train = pd.concat([X_train, X_test.loc[list(mutation_pair)]])
            y_train = pd.concat([y_train, y_test.loc[list(mutation_pair)]])

            # Remove from test set
            X_test = X_test.drop(mutation_pair)
            y_test = y_test.drop(mutation_pair)

            pairs_leaked += 1  # Increment the count of leaked pairs

            # Retrain and evaluate
            model.fit(X_train, y_train)
            y_pred = model.predict(X_validation)
            pcc_post = np.corrcoef(y_validation, y_pred)[0, 1]
            mae_post = np.mean(np.abs(y_pred - y_validation))

            # Track leakage performance
            leakage_tracking[pdb_id]['pcc'].append(pcc_post)
            leakage_tracking[pdb_id]['mae'].append(mae_post)
            leakage_tracking[pdb_id]['leaked_pairs'].append(pairs_leaked)
            print(f"{pdb_id} PCC after leaking {pairs_leaked} pairs: {pcc_post}")

        pccs[pdb_id + "_postleak"] = pcc_post
        print(f"{pdb_id} Post-leak PCC: {pcc_post}")

        pcc_improvements[pdb_id] = pcc_post - pcc_pre
        print(f"{pdb_id} PCC Improvement: {pcc_improvements[pdb_id]}")

        feature_importances[pdb_id] = model.feature_importances_

        if pcc_improvements[pdb_id] > 0:
            improved_pdbs.append(pdb_id)
        else:
            worse_pdbs.append(pdb_id)

        if i == 21:
            break

    return leakage_tracking, pccs

def plot_pca_complex(X, complex_ids, n_components=2):
    # Initialize PCA
    pca = PCA(n_components=n_components)

    # Fit PCA on all data  
    X_reduced = pca.fit_transform(X)
    X_reduced_df = pd.DataFrame(X_reduced, index=X.index)

    # Create the plot
    plt.figure(figsize=(12, 6))

    # Plot all other complexes in dark grey
    other_complexes_indices = ~X_reduced_df.index.str.startswith(tuple(complex_ids))
    other_data = X_reduced_df[other_complexes_indices]
    plt.scatter(other_data[0], other_data[1], label='Other Complexes', color='darkgrey', alpha=0.5)

    # Plot each specified complex with a different color
    colors = plt.colormaps.get_cmap('viridis')(np.linspace(0, 1, len(complex_ids))) 
    for i, complex_id in enumerate(complex_ids):
        complex_indices = X_reduced_df.index.str.startswith(complex_id)
        complex_data = X_reduced_df[complex_indices]
        plt.scatter(complex_data[0], complex_data[1], label=complex_id, color=colors[i], alpha=0.7)

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA of Feature Space (Highlighted Complexes)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)  # Place legend outside the plot
    plt.grid(True)
    plt.tight_layout()

def plot_umap_complex(X, complex_ids, n_components=2, random_state=None):

    # Initialize UMAP
    reducer = umap.UMAP(n_components=n_components, random_state=random_state)

    # Fit UMAP on all data
    X_reduced = reducer.fit_transform(X)
    X_reduced_df = pd.DataFrame(X_reduced, index=X.index)

    # Create the plot
    plt.figure(figsize=(10, 8))

    # Plot all other complexes in dark grey
    other_complexes_indices = ~X_reduced_df.index.str.startswith(tuple(complex_ids))
    other_data = X_reduced_df[other_complexes_indices]
    plt.scatter(other_data[0], other_data[1], label='Other Complexes', color='darkgrey', alpha=0.5)

    # Plot each specified complex with a different color
    colors = plt.colormaps.get_cmap('viridis')(np.linspace(0, 1, len(complex_ids)))  # Generate a colormap
    for i, complex_id in enumerate(complex_ids):
        complex_indices = X_reduced_df.index.str.startswith(complex_id)
        complex_data = X_reduced_df[complex_indices]
        plt.scatter(complex_data[0], complex_data[1], label=complex_id, color=colors[i], alpha=0.7)

    plt.xlabel('UMAP Component 1')
    plt.ylabel('UMAP Component 2')
    plt.title('UMAP of Feature Space (Highlighted Complexes)')
    plt.legend()
    plt.grid(True)

def plot_tsne_complex(X, complex_ids, n_components=2, random_state=None):

    # Initialize t-SNE
    tsne = TSNE(n_components=n_components, random_state=random_state)

    # Fit t-SNE on all data
    X_reduced = tsne.fit_transform(X)
    X_reduced_df = pd.DataFrame(X_reduced, index=X.index)

    # Create the plot
    plt.figure(figsize=(11, 6))

    # Plot all other complexes in dark grey
    other_complexes_indices = ~X_reduced_df.index.str.startswith(tuple(complex_ids))
    other_data = X_reduced_df[other_complexes_indices]
    plt.scatter(other_data[0], other_data[1], label='Other Complexes', color='darkgrey', alpha=0.5)

    # Plot each specified complex with a different color
    # colors = plt.colormaps.get_cmap('viridis')(np.linspace(0, 1, len(complex_ids)))  # Generate a colormap
    colors = plt.colormaps.get_cmap('rainbow')(np.linspace(0, 1, len(complex_ids)))  # Generate a colormap
    for i, complex_id in enumerate(complex_ids):
        complex_indices = X_reduced_df.index.str.startswith(complex_id)
        complex_data = X_reduced_df[complex_indices]
        plt.scatter(complex_data[0], complex_data[1], label=complex_id, color=colors[i], alpha=0.7)

    plt.xlabel('t-SNE Component 1')
    plt.ylabel('t-SNE Component 2')
    plt.title('t-SNE of Feature Space (Highlighted Complexes)')
    plt.legend()
    plt.grid(True)

if __name__ == '__main__':
    RANDOM = 167

    # load the features
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
    # print the shape of the delta features
    print('all delta features total', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # permute mutations
    delta_features = permute_mutations(delta_features)
    print('all delta features permuted', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # reverse mutations
    delta_features = reverse_mutations(delta_features)
    print('all delta features permuted and reversed', len(delta_features), ', ', len(delta_features[list(delta_features.keys())[0]]))

    # convert to dataframe
    delta_features_df = pd.DataFrame(delta_features).T
    # print(delta_features_df.shape)
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

    # keep only select features from recursive feature elimination
    delta_features_df = keep_selected_features(delta_features_df, '../data/selected_features.txt')
    print(delta_features_df.shape, " after keeping select features")

    # resample overrepresented complexes to balance the dataset
    delta_features_df = balance_sample(delta_features_df, threshold=b)
    print(delta_features_df.shape, " after balancing the data")

    # # write the balanced data to a file
    # delta_features_df.to_csv('selected_delta_features.csv')
    # delta_features_df = pd.read_csv('selected_delta_features.csv', index_col=0)

    X = delta_features_df.drop(columns=['ddG'])
    y = delta_features_df['ddG']

    # train and test the model
    leakage_tracking, pccs = test_leakage(X, y, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, max_features=max_features)
    
    # keep only the good pdbs
    # # good pdbs from the single state model
    # good_pdbs = ['1DAN', '4P5T', '3BT1', '4P23', '1LFD', '1IAR', '1A22', '1KTZ', '1CHO', '3QIB']
    # good pdbs from the ensemble model
    good_pdbs = ['1DAN', '4P5T', '3BT1', '4P23', '1LFD', '1IAR', '1A22', '3S9D', '1BRS', '1EMV']
    bad_pdbs = [pdb for pdb in leakage_tracking if pdb not in good_pdbs]
    # remove 36C0
    # bad_pdbs.remove('3C60')
    bad_pdbs.remove('2AK4')
    bad_pdbs.remove('1JTG')
    print("Good PDBs:", good_pdbs)
    print("Bad PDBs:", bad_pdbs)

    # filter out the bad pdbs from the leakage tracking
    leakage_tracking_bad = {k: v for k, v in leakage_tracking.items() if k in bad_pdbs}
    leakage_tracking = {k: v for k, v in leakage_tracking.items() if k in good_pdbs}

    # Prepare and save good PDBs leakage tracking
    good_pcc_df = pd.DataFrame({k: v['pcc'] for k, v in leakage_tracking.items()}).T
    good_mae_df = pd.DataFrame({k: v['mae'] for k, v in leakage_tracking.items()}).T
    good_pcc_df.columns = [f'pcc_{i}' for i in range(good_pcc_df.shape[1])]
    good_mae_df.columns = [f'mae_{i}' for i in range(good_mae_df.shape[1])]
    good_leakage_df = pd.concat([good_pcc_df, good_mae_df], axis=1)
    good_leakage_df.to_csv('leakage_tracking_good.csv')

    # Prepare and save bad PDBs leakage tracking
    bad_pcc_df = pd.DataFrame({k: v['pcc'] for k, v in leakage_tracking_bad.items()}).T
    bad_mae_df = pd.DataFrame({k: v['mae'] for k, v in leakage_tracking_bad.items()}).T
    bad_pcc_df.columns = [f'pcc_{i}' for i in range(bad_pcc_df.shape[1])]
    bad_mae_df.columns = [f'mae_{i}' for i in range(bad_mae_df.shape[1])]
    bad_leakage_df = pd.concat([bad_pcc_df, bad_mae_df], axis=1)
    bad_leakage_df.to_csv('leakage_tracking_bad.csv')

    # Pre- and post-leak PCC distribution plot
    fig, axs = plt.subplots(2, 1, figsize=(8,10))
    preleak_pccs = [pccs[k] for k in pccs if k.endswith("_preleak")]
    postleak_pccs = [pccs[k] for k in pccs if k.endswith("_postleak")]

    # Incremental leakage performance plots (1 mutation pair at a time)
    # axs[0].figure(figsize=(10, 6))  # Adjust figure size for better legend placement
    num_complexes = len(leakage_tracking)
    # colors = plt.colormaps.get_cmap('viridis')(np.linspace(0, 1, num_complexes)) 
    # make colors the rainbow colors
    colors = plt.colormaps.get_cmap('rainbow')(np.linspace(0, 1, num_complexes))

    for i, (pdb_id, data) in enumerate(leakage_tracking.items()):
        axs[0].plot(data['leaked_pairs'], data['pcc'], label=pdb_id, color=colors[i])

    axs[0].set_xlabel("Number of Mutation Pairs Leaked")
    axs[0].set_ylabel("PCC")
    axs[0].set_title("Incremental PCC Change per Complex")
    # set ylimits to be 0-1
    axs[0].set_ylim(0.3, 0.85)
    # axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)  # Place legend outside the plot
    # fig.tight_layout()  # Adjust layout to make room for the legend

    # plt.figure(figsize=(10, 6))  # Adjust figure size for better legend placement
    for i, (pdb_id, data) in enumerate(leakage_tracking.items()):
        axs[1].plot(data['leaked_pairs'], data['mae'], label=pdb_id, color=colors[i])

    axs[1].set_xlabel("Number of Mutation Pairs Leaked")
    axs[1].set_ylabel("MAE")
    axs[1].set_title("Incremental MAE Change per Complex")
    axs[1].set_ylim(0.3, 2.25)

    fig.tight_layout()
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1, 0.5), ncol=1)

    # subplot adjust
    plt.subplots_adjust(left=0.1, right=0.75, top=0.95, bottom=0.1, wspace=0.2, hspace=0.32)


    # plot the bad pdbs leakage tracking
    fig, axs = plt.subplots(2, 1, figsize=(8,10))
    num_complexes = len(leakage_tracking_bad)
    # colors = plt.colormaps.get_cmap('viridis')(np.linspace(0, 1, num_complexes))
    colors = plt.colormaps.get_cmap('rainbow')(np.linspace(0, 1, num_complexes))

    for i, (pdb_id, data) in enumerate(leakage_tracking_bad.items()):
        axs[0].plot(data['leaked_pairs'], data['pcc'], label=pdb_id, color=colors[i])

    axs[0].set_xlabel("Number of Mutation Pairs Leaked")
    axs[0].set_ylabel("PCC")
    axs[0].set_title("Incremental PCC Change per Complex")
    axs[0].set_ylim(0.3, 0.85)
    # axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)  # Place legend outside the plot
    # fig.tight_layout()  # Adjust layout to make room for the legend

    # plt.figure(figsize=(10, 6))  # Adjust figure size for better legend placement
    for i, (pdb_id, data) in enumerate(leakage_tracking_bad.items()):
        axs[1].plot(data['leaked_pairs'], data['mae'], label=pdb_id, color=colors[i])
        
    axs[1].set_xlabel("Number of Mutation Pairs Leaked")
    axs[1].set_ylabel("MAE")
    axs[1].set_title("Incremental MAE Change per Complex")
    axs[1].set_ylim(0.3, 2.25)
    
    # fig.tight_layout()
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1, 0.5), ncol=1)

    # subplot adjust
    plt.subplots_adjust(left=0.1, right=0.75, top=0.95, bottom=0.1, wspace=0.2, hspace=0.32)


    # plot pca
    # plot_pca_complex(X, good_pdbs)
    # plot_tsne_complex(X, good_pdbs)

    # plot pca
    # plot_pca_complex(X, bad_pdbs)
    # plot_tsne_complex(X, bad_pdbs)

    # Calculate average % increase of PCCs in good pdbs
    pcc_increase_good = []
    for pdb_id in good_pdbs:
        preleak_pcc = pccs[pdb_id + "_preleak"]
        postleak_pcc = pccs[pdb_id + "_postleak"]
        pcc_increase = ((postleak_pcc - preleak_pcc) / abs(preleak_pcc)) * 100
        pcc_increase_good.append(pcc_increase)
    avg_pcc_increase_good = np.mean(pcc_increase_good)

    # Calculate average % increase of PCCs in bad pdbs
    pcc_increase_bad = []
    for pdb_id in bad_pdbs:
        preleak_pcc = pccs[pdb_id + "_preleak"]
        postleak_pcc = pccs[pdb_id + "_postleak"]
        pcc_increase = ((postleak_pcc - preleak_pcc) / abs(preleak_pcc)) * 100
        pcc_increase_bad.append(pcc_increase)
    avg_pcc_increase_bad = np.mean(pcc_increase_bad)

    # Calculate average % decrease in MAE for good pdbs
    mae_decrease_good = []
    for pdb_id, data in leakage_tracking.items():
        preleak_mae = data['mae'][0]
        postleak_mae = data['mae'][-1]
        mae_decrease = ((preleak_mae - postleak_mae) / abs(preleak_mae)) * 100
        mae_decrease_good.append(mae_decrease)
    avg_mae_decrease_good = np.mean(mae_decrease_good)

    # Calculate average % decrease in MAE for bad pdbs
    mae_decrease_bad = []
    for pdb_id, data in leakage_tracking_bad.items():
        preleak_mae = data['mae'][0]
        postleak_mae = data['mae'][-1]
        mae_decrease = ((preleak_mae - postleak_mae) / abs(preleak_mae)) * 100
        mae_decrease_bad.append(mae_decrease)
    avg_mae_decrease_bad = np.mean(mae_decrease_bad)

    print(f"Average % increase of PCC in good pdbs: {avg_pcc_increase_good:.2f}%")
    print(f"Average % increase of PCC in bad pdbs: {avg_pcc_increase_bad:.2f}%")
    print(f"Average % decrease in MAE for good pdbs: {avg_mae_decrease_good:.2f}%")
    print(f"Average % decrease in MAE for bad pdbs: {avg_mae_decrease_bad:.2f}%")

    # # Feature Space Expansion Analysis
    # X = delta_features_df.drop(columns=['ddG'])
    
    # # Scale the features
    # scaler = StandardScaler()
    # X_scaled = scaler.fit_transform(X)
    # X_scaled = pd.DataFrame(X_scaled, index=X.index, columns=X.columns)

    # # Calculate the center of the feature space using all training data
    # feature_space_center = X_scaled.mean(axis=0)

    # def calculate_expansion(pdb_ids, X_scaled, feature_space_center):
    #     distances = {}
    #     for pdb_id in pdb_ids:
    #         # Get the indices for the current PDB ID
    #         pdb_indices = X_scaled.index.str.startswith(pdb_id)
            
    #         # Calculate Euclidean distances from the center for the current PDB ID
    #         pdb_distances = euclidean_distances(X_scaled[pdb_indices], feature_space_center.values.reshape(1, -1))
    #         distances[pdb_id] = np.mean(pdb_distances), np.std(pdb_distances)
    #     return distances

    # # Calculate expansion for good PDBs
    # good_distances = calculate_expansion(good_pdbs, X_scaled, feature_space_center)
    # mean_distance_good = np.mean([good_distances[pdb][0] for pdb in good_distances])
    # std_distance_good = np.std([good_distances[pdb][0] for pdb in good_distances])
    # print(f"Good PDBs - Avg Distance from Center: {mean_distance_good:.2f}, Std Dev: {std_distance_good:.2f}")

    # # Calculate expansion for bad PDBs
    # bad_distances = calculate_expansion(bad_pdbs, X_scaled, feature_space_center)
    # mean_distance_bad = np.mean([bad_distances[pdb][0] for pdb in bad_distances])
    # std_distance_bad = np.std([bad_distances[pdb][0] for pdb in bad_distances])
    # print(f"Bad PDBs - Avg Distance from Center: {mean_distance_bad:.2f}, Std Dev: {std_distance_bad:.2f}")

    # # sort the good and bad distances
    # good_distances = dict(sorted(good_distances.items(), key=lambda item: item[1][0], reverse=True))
    # bad_distances = dict(sorted(bad_distances.items(), key=lambda item: item[1][0], reverse=True))

    # # Prepare data for bar plot
    # pdb_ids = list(good_distances.keys()) + list(bad_distances.keys())
    # distances = [good_distances[pdb][0] for pdb in good_distances] + [bad_distances[pdb][0] for pdb in bad_distances]
    # colors = ['blue'] * len(good_distances) + ['red'] * len(bad_distances)

    # # Create bar plot
    # plt.figure(figsize=(12, 6))
    # plt.bar(pdb_ids, distances, color=colors)
    # plt.xlabel("PDB ID")
    # plt.ylabel("Average Distance from Feature Space Center")
    # plt.title("Average Distance from Feature Space Center for Good and Bad PDBs")
    # plt.xticks(rotation=45, ha="right")
    # plt.tight_layout()

    plt.show()