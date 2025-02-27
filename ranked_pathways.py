import os
import argparse
import pandas as pd
import anndata as ad

def load_data(h5ad_file):
    return ad.read_h5ad(h5ad_file)

def filter_data(adata, ligrec_file):
    df = pd.read_csv(ligrec_file, sep="\t", header=None, names=["ligand", "receptor", "pathway"])
    ligands = df["ligand"].tolist()
    receptors = df["receptor"].tolist()
    pathways = df["pathway"].tolist()
    
    combined_strings = [f"commot-user_database-{ligand}-{receptor}" for ligand, receptor in zip(ligands, receptors)]
    
    mean_values = []
    for key in combined_strings:
        if key in adata.obsp:
            sparse_matrix = adata.obsp[key]
            mean_value = sparse_matrix.mean()
            mean_values.append(mean_value)
        else:
            mean_values.append(0)  
    
    df_mean = pd.DataFrame({
        "pathway": pathways,
        "mean_value": mean_values
    })
    grouped_means = df_mean.groupby("pathway").mean().reset_index()
    return grouped_means, ligands, receptors, pathways, mean_values

def main(args):
    adata = load_data(args.adata_file)
    
    grouped_means, _, _, _, _ = filter_data(adata, args.ligrec_file)
    
    sorted_means = grouped_means.sort_values(by='mean_value', ascending=False)
    
    sorted_means.to_csv(args.output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print top 20 ligand-receptor-pathway combinations with the highest mean values.')
    parser.add_argument('--adata_file', type=str, required=True, help='Path to the AnnData file.')
    parser.add_argument('--ligrec_file', type=str, required=True, help='Path to the ligand-receptor pairs file.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the output file.')
    args = parser.parse_args()
    main(args)

