import pandas as pd
from matplotlib import pyplot as plt

# Load the data
data = pd.read_csv('correlations.csv')

hyperparameters = ['N', 'B', 'H', 'N_ESTIMATORS', 'MAX_DEPTH', 'MIN_SAMPLES_SPLIT', 'MIN_SAMPLES_LEAF', 'MAX_FEATURES']
objective = ["PCC"]

# Convert MAX_FEATURES to numerical if it's not already
if data['MAX_FEATURES'].dtype == 'object':
    data['MAX_FEATURES'] = data['MAX_FEATURES'].astype('category').cat.codes

# Calculate mean PCC for each hyperparameter value
mean_pcc_per_param = {}
for param in hyperparameters:
    mean_pcc = data.groupby(param)['PCC'].mean().sort_values(ascending=False)
    mean_pcc_per_param[param] = mean_pcc
    print(f"Mean PCC for each value of {param}:\n{mean_pcc}\n")

# Identify the best combination of parameters
best_params = data.sort_values(by='PCC', ascending=False).iloc[0]

print("\nBest Combination of Hyperparameters:")
print(best_params)

# make difference heatmaps for each pair of hyperparameters
for i in range(len(hyperparameters)):
    for j in range(i+1, len(hyperparameters)):
        # make a heatmap of the difference in PCC between the two hyperparameters
        fig, ax = plt.subplots()
        ax.set_title(f"{hyperparameters[i]} vs {hyperparameters[j]}")
        ax.set_xlabel(hyperparameters[i])
        ax.set_ylabel(hyperparameters[j])
        
        # Pivot the data, handling duplicates by averaging
        pivot_table = data.pivot_table(index=hyperparameters[j], columns=hyperparameters[i], values=objective, aggfunc='mean')
        
        heatmap = ax.imshow(pivot_table, cmap='coolwarm')
        plt.colorbar(heatmap)
        plt.savefig(f"{hyperparameters[i]}_vs_{hyperparameters[j]}.png")
        plt.close()

# plt.show()