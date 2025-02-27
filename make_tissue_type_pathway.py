import anndata as ad
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="Process AnnData objects and filter spatial coordinates.")
parser.add_argument("--adata1_path", required=True, help="Path to the first AnnData file.")
parser.add_argument("--adata2_path", required=True, help="Path to the second AnnData file.")
parser.add_argument("--csv_path", required=True, help="Path to the CSV file.")
parser.add_argument("--pathway_file", required=True, help="Path to the file containing pathway names.")
parser.add_argument("--output_dir", required=True, help="Directory to save the output CSV files.")
parser.add_argument("--tissue_type", required=True, help="Tissue type to include in the output file name.")
args = parser.parse_args()

adata1_path = args.adata1_path
adata2_path = args.adata2_path
csv_path = args.csv_path
pathway_file = args.pathway_file
output_dir = args.output_dir
tissue_type = args.tissue_type

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

adata1 = ad.read_h5ad(adata1_path)
adata2 = ad.read_h5ad(adata2_path)
data = pd.read_csv(csv_path)

with open(pathway_file, 'r') as f:
    pathways = [line.strip() for line in f]

spatial_coords = adata1.obsm["spatial"]

valid_obs_names = set(data.iloc[:, 0])

for pathway in pathways:
    obsp_key = f"commot-user_database-{pathway}"
    print(f"Processing {obsp_key}...")
    
    if obsp_key not in adata2.obsp:
        print(f"{obsp_key} not found in adata2.obsp. Skipping.")
        continue
    
    matrix = adata2.obsp[obsp_key]
    nonzero_indices = list(zip(*matrix.nonzero()))
    
    filtered_indices = [
        (sender_idx, receiver_idx)
        for sender_idx, receiver_idx in nonzero_indices
        if adata1.obs_names[sender_idx] in valid_obs_names and adata1.obs_names[receiver_idx] in valid_obs_names
    ]
    
    if filtered_indices:
        results = []
        for sender_idx, receiver_idx in filtered_indices:
            results.append({
                "Sender": adata1.obs_names[sender_idx],
                "Receiver": adata1.obs_names[receiver_idx],
                "Sender_x": spatial_coords[sender_idx][0],
                "Sender_y": spatial_coords[sender_idx][1],
                "Receiver_x": spatial_coords[receiver_idx][0],
                "Receiver_y": spatial_coords[receiver_idx][1],
            })
        
        results_df = pd.DataFrame(results)
        
        coords_map = data.set_index(data.columns[0])[['xcoord', 'ycoord']].to_dict('index')
        
        results_df['Sender_x'] = results_df['Sender'].map(lambda x: coords_map.get(x, {}).get('xcoord'))
        results_df['Sender_y'] = results_df['Sender'].map(lambda x: coords_map.get(x, {}).get('ycoord'))
        results_df['Receiver_x'] = results_df['Receiver'].map(lambda x: coords_map.get(x, {}).get('xcoord'))
        results_df['Receiver_y'] = results_df['Receiver'].map(lambda x: coords_map.get(x, {}).get('ycoord'))
        
        results_df.drop_duplicates(inplace=True)
        
        output_file = os.path.join(output_dir, f"{tissue_type}_{pathway}.csv")
        results_df.to_csv(output_file, index=False)

        print(f"Filtered coordinates for {obsp_key} saved to {output_file}")
    else:
        print(f"No valid sender-receiver pairs found for {obsp_key}. Skipping.")

