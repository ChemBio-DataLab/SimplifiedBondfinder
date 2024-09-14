# 
# BSD 3-Clause License
# Copyright (c) 2023-2024 Sophia Bazzi and Sharareh Sayyad
# All rights reserved.

import pandas as pd
from umap import UMAP
from sklearn.cluster import MeanShift
from sklearn.metrics import silhouette_score, davies_bouldin_score
from itertools import combinations
import numpy as np
import hashlib
from sklearn.preprocessing import MinMaxScaler

# Your list of reference pdb_filenames
pdb_filenames = [
"6cb0_A_CYS_SG_257_A_LYS_NZ_255.pdb", "6cb0_B_CYS_SG_257_B_LYS_NZ_255.pdb", "2ccm_A_CYS_SG_24_A_LYS_NZ_41.pdb", "2ccm_B_CYS_SG_24_B_LYS_NZ_41.pdb" #NOS from Nature paper
]

# Define the hash_set function for consistent hashing of sets
def hash_set(set_data):
    try:
        # Check if set_data is already a set or iterable
        if isinstance(set_data, (list, tuple, set)):
            # Convert to set if not already a set
            set_values = set(set_data)
        elif isinstance(set_data, str):
            # Parse string representation of set into a set of values
            set_values = set(eval(set_data))  # Use eval with caution to parse literal expressions
        else:
            # If set_data is a single value (e.g., int), convert it to a single-element set
            set_values = {set_data}

        # Sort the set values and convert to tuple for hashing
        sorted_tuple = tuple(sorted(set_values))
        hashed_value = hash(sorted_tuple)

        return hashed_value

    except Exception as e:
        print(f"Error processing set_data: {e}")
        return None

# Load your data
data = pd.read_csv('CYS_non_N_CYS_LYS_No_back.txt', sep='\t')

# Convert 'ngb_Res_toS' and 'ngb_Res_toN' to consistent hash values in the original DataFrame
if 'ngb_Res_toS' in data.columns:
    data['ngb_Res_toS'] = data['ngb_Res_toS'].apply(hash_set)
    # Apply z-score normalization to 'ngb_Res_toS'
    data['ngb_Res_toS_normalized'] = (data['ngb_Res_toS'] - data['ngb_Res_toS'].mean()) / data['ngb_Res_toS'].std()

if 'ngb_Res_toN' in data.columns:
    data['ngb_Res_toN'] = data['ngb_Res_toN'].apply(hash_set)
    # Apply z-score normalization to 'ngb_Res_toN'
    data['ngb_Res_toN_normalized'] = (data['ngb_Res_toN'] - data['ngb_Res_toN'].mean()) / data['ngb_Res_toN'].std()

# Save the modified dataset with hash columns to a new file (tab-separated with updated header)
output_file_path = 'CYS_non_N_CYS_LYS_No_back_hashes.csv'
#data.to_csv(output_file_path, sep='\t', index=False)
data.to_csv(output_file_path, sep='\t', index=False)
print(f"Modified data saved to {output_file_path}")

# Save 'filename' column for later use
filenames = data['filename']  # Save the 'filename' column

columns_to_drop = [
    'Chain_S', 'Res_S', 'S_id', 'CS_id', 'Res_S_num', 'Chain_N', 'Res_N', 'N_id', 'CN_id', 'Res_N_num',
    'greenBlob', 'Res_S_num', 'Close_contact', 'Res_CS', 'bfac_S', 'bfac_N',
    'distanceCS', 'ngb_Res_toN', 'ngb_Res_toS',
    'bfac_CN', 'bfac_CS',
    'occ_S', 'occ_N', 'occu_CN', 'occu_CS',
]
data = data.drop(columns=columns_to_drop)
# Define the feature_space excluding 'filename'
feature_space = [col for col in data.columns if col != 'filename']

# Write the results to files
with open("clusters_with_common_filenames.txt", "w") as detailed_file, \
        open("clusters_summary.txt", "w") as summary_file:
    for num_features in range(1, len(feature_space) + 1):  # Iterate over feature_space length
        for feature_combo in combinations(feature_space, num_features):
            subset_data = data[list(feature_combo)].copy()  # Subset data with selected features

            if len(subset_data.columns) < 3:
                continue  # Skip if subset_data contains fewer than 2 features

            try:
                # Set the values for min_dist and n_neighbors
                print(feature_space)
                min_dist = 0.0
                n_neighbors = 7

                # Apply UMAP for dimensionality reduction
                umap_model = UMAP(n_components=3, random_state=42, min_dist=min_dist, n_neighbors=n_neighbors)
                umap_result = umap_model.fit_transform(subset_data)

                # Perform clustering using MeanShift
                np.random.seed(42)
                mean_shift = MeanShift()
                cluster_labels = mean_shift.fit_predict(umap_result)

                # Check if there's more than one unique label
                if len(set(cluster_labels)) <= 1:
                    continue  # Skip this combination if only one cluster is found

                # Calculate silhouette score
                silhouette_avg = silhouette_score(umap_result, cluster_labels)

                # Calculate Davies–Bouldin index
                dbi = davies_bouldin_score(umap_result, cluster_labels)

            except Exception as e:
                print("An error occurred:", e)
                continue  # Continue with the next iteration of the loop

            # Write results to detailed file
            detailed_file.write(f"Features: {', '.join(feature_combo)}\n")
            detailed_file.write(f"Silhouette Score: {silhouette_avg}\n")
            detailed_file.write(f"Davies–Bouldin Index: {dbi}\n")
            detailed_file.write('-' * 50 + '\n')

            # Print filenames that belong to each cluster
            cluster_counts = {}
            for cluster in set(cluster_labels):
                cluster_indices = [i for i, label in enumerate(cluster_labels) if label == cluster]
                cluster_file_names = filenames.iloc[cluster_indices]  # Get filenames for this cluster
                cluster_file_names_processed = list(cluster_file_names)
                common_filenames = set(cluster_file_names_processed) & set(pdb_filenames)  # Common with PDB filenames

                detailed_file.write(f"Cluster {cluster} File Names:\n")
                detailed_file.write('\n'.join([f"Cluster {cluster}: {filename}" for filename in cluster_file_names_processed]) + '\n\n')
                detailed_file.write("Common with PDB Filenames:\n")
                detailed_file.write('\n'.join([f"Cluster {cluster}: {common_filename}" for common_filename in common_filenames]) + '\n\n')
                detailed_file.write(f"Count of Common with PDB Filenames for cluster {cluster}: {len(common_filenames)}\n")
                detailed_file.write('-' * 100 + '\n')

                # Update cluster count
                cluster_counts[cluster] = len(cluster_file_names_processed)

            # Write results to summary file
            summary_file.write(f"Features: {', '.join(feature_combo)}\n")
            summary_file.write(f"Silhouette Score: {silhouette_avg}\n")
            summary_file.write(f"Davies–Bouldin Index: {dbi}\n")
            summary_file.write('-' * 50 + '\n')

            for cluster, count in cluster_counts.items():
                cluster_indices = [i for i, label in enumerate(cluster_labels) if label == cluster]
                cluster_file_names = filenames.iloc[cluster_indices]  # Get filenames for this cluster
                common_filenames = set(cluster_file_names) & set(pdb_filenames)  # Common with PDB filenames

                summary_file.write(f"Cluster {cluster}: {count} filenames\n")
                summary_file.write(f"Count of Common with PDB Filenames for cluster {cluster}: {len(common_filenames)}\n")
                summary_file.write('-' * 100 + '\n')

