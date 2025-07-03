import pandas as pd

def split_data(file):
    # split the features into original and mutated features
    df = pd.read_csv(file, index_col=0)
    # get the original features
    original_features = df[df.index.str.contains('_filled')]
    # get the mutated features
    mutated_features = df[~df.index.str.contains('_filled')]
    # write the original features to a file
    original_features.to_csv('wildtype_features.csv')
    # write the mutated features to a file
    mutated_features.to_csv('single_state_eppi_features.csv')

def get_overlapping_mutations(file1, file2):
    # get the overlapping mutations between two files
    df1 = pd.read_csv(file1, index_col=0)
    df2 = pd.read_csv(file2, index_col=0)
    # get the overlapping mutations
    overlapping_mutations = df1.index.intersection(df2.index)
    # write the overlapping mutations to two subset files
    df1.loc[overlapping_mutations].to_csv(file1.split('.')[0] + '_subset.csv')
    df2.loc[overlapping_mutations].to_csv(file2.split('.')[0] + '_subset.csv')
    return overlapping_mutations


if __name__ == '__main__':
    # split_data('all_data.csv')
    get_overlapping_mutations('wildtype_features.csv', 'ensemble_eppi_features.csv')
