import os
import argparse
import pandas as pd
import anndata as ad
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
import numpy as np
import plotly.graph_objects as go

def load_data(h5ad_file):
    """Load data from an AnnData file."""
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
        "ligand": ligands,
        "receptor": receptors,
        "pathway": pathways,
        "mean": mean_values
    })
     
    df_mean_filtered = df_mean[df_mean["mean"] != 0]
    
    return df_mean_filtered

def rank_ligand_receptor(df_mean):
    top_ligand_receptor = df_mean.nlargest(10, 'mean')
    
    top_ligand_receptor['z_score'] = zscore(top_ligand_receptor['mean'])
    
    min_z_score = top_ligand_receptor['z_score'].min()
    max_z_score = top_ligand_receptor['z_score'].max()
    top_ligand_receptor['normalized_score'] = ((top_ligand_receptor['z_score'] - min_z_score) / (max_z_score - min_z_score)) * 99 + 1

    return top_ligand_receptor

def create_sankey_plot(top_ligand_receptor, output_file):
    pathway_scores = top_ligand_receptor.groupby('pathway')['mean'].mean().reset_index()
    pathway_scores['z_score'] = zscore(pathway_scores['mean'])
    min_z_score = pathway_scores['z_score'].min()
    max_z_score = pathway_scores['z_score'].max()
    pathway_scores['normalized_score'] = ((pathway_scores['z_score'] - min_z_score) / (max_z_score - min_z_score)) * 99 + 1
    print(pathway_scores)
    sankey_df = 

def create_sankey_plot1(top_ligand_receptor, output_file):
    pathway_scores = top_10.groupby('pathway')['normalized_score'].mean()
    min_z_score = zscore(pathway_scores).iloc[0]
    max_z_score = zscore(pathway_scores).iloc[-1]
    new_normalized_scores = ((zscore(pathway_scores) - min_z_score) / (max_z_score - min_z_score)) * 99 + 1
    pathway_score_mapping = dict(zip(pathway_scores.index, new_normalized_scores))
    
    sankey_df = pd.DataFrame(columns=['source', 'target', 'value'])
    
    sources = pd.concat([top_10['ligand'], top_10['receptor']]).unique()
    
    targets = pd.concat([top_10['receptor'], top_10['pathway']]).unique()
    print(top_10)
    print(sources)
    rows = []
    
    for ligand, receptor, pathway, normalized_score in top_10[['ligand', 'receptor', 'pathway', 'normalized_score']].itertuples(index=False):
        for target in targets:
            if target in top_10['receptor'].unique():
                value = normalized_score
            else:
                value = pathway_score_mapping[target]
            rows.append({'source': receptor, 'target': pathway, 'value': value})
    
    sankey_df = pd.concat([pd.DataFrame([row]) for row in rows], ignore_index=True)
    
    print(sankey_df)
    
    labels = list(sources) + list(set(targets) - set(sources))
    
    label_indices = {label: idx for idx, label in enumerate(labels)}
    
    sankey_data = {
        'source': [label_indices[row['source']] for row in rows],
        'target': [label_indices[row['target']] for row in rows],
        'value': [row['value'] for row in rows]
    }
    print(labels)
    print(sankey_data)
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color="blue"),
        link=dict(
            source=sankey_data['source'],
            target=sankey_data['target'],
            value=sankey_data['value']))])

def main(args):
    adata = load_data(args.adata_file)
    
    df_mean = filter_data(adata, args.ligrec_file)
    
    output_file = "ligand_receptor_pathway.tsv"
    df_mean.to_csv(output_file, sep='\t', index=False)
    print(f"Ligand-receptor-pathway-mean has been saved to {output_file}")
    
    top_ligand_receptor = rank_ligand_receptor(df_mean)
    
    output_file = "top_ligand_receptor_pathway.tsv"
    top_ligand_receptor.to_csv(output_file, sep='\t', index=False)
    print(f"Ligand-receptor-pathway-mean-z_score-normalized_score has been saved to {output_file}")
    
    create_sankey_plot(top_ligand_receptor, "sankey_plot.html")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print top 20 ligand-receptor-pathway combinations with the highest mean values.')
    parser.add_argument('--adata_file', type=str, required=True, help='Path to the AnnData file.')
    parser.add_argument('--ligrec_file', type=str, required=True, help='Path to the ligand-receptor pairs file.')
    args = parser.parse_args()
    main(args)

