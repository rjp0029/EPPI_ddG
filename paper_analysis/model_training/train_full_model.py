import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import joblib

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

    # get X and y
    y = delta_features_df['ddG']
    X = delta_features_df.drop(columns=['ddG'])

    # train model
    model = RandomForestRegressor(random_state=RANDOM, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, max_features=max_features)
    model.fit(X, y)

    # save the model
    joblib.dump(model, 'ensemble_full.pkl')
