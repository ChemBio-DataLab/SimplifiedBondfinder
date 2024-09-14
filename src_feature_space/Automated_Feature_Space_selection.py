# 
# BSD 3-Clause License
# Copyright (c) 2023-2024 Sophia Bazzi and Sharareh Sayyad
# All rights reserved.

def extract_features(file_path, output_file):
    features = []
    current_feature = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("Features:"):
                current_feature = {'name': line.split(": ")[1], 'silhouette_score': None, 'clusters': []}
            elif line.startswith("Silhouette Score:"):
                current_feature['silhouette_score'] = float(line.split(": ")[1])
            elif line.startswith("Cluster"):
                cluster_info = line.strip() + " " + lines[i+1].strip()
                current_feature['clusters'].append(cluster_info)
            elif line.startswith("Count of Common with PDB Filenames"):
                count_common = line.split(": ")[1]
                current_feature['clusters'][-1] += ", " + count_common
                if current_feature['silhouette_score'] is not None and current_feature['silhouette_score'] >= 0.6 and int(count_common) == 81:
                    features.append(current_feature)
            i += 1
    
    # Write the selected features to the output file
    with open(output_file, 'w') as output:
        for feature in features:
            output.write("Feature: {}\n".format(feature['name']))
            output.write("Silhouette Score: {}\n".format(feature['silhouette_score']))
            output.write("Clusters:\n")
            for cluster in feature['clusters']:
                output.write(cluster + "\n")
            output.write("\n")

file_path = "clusters_summary.txt"
output_file = "selected_features.txt"
extract_features(file_path, output_file)
print("Selected features have been written to", output_file)

