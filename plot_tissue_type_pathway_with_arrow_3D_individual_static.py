import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_cells_with_connections_static(
    all_coords,
    highlight_coords,
    sender_receiver_data,
    tissue_type_name,
    pathway_name,
    static_file=None,
    width=10,
    height=10,
    alpha_all=0.8,
    alpha_highlight=0.9,
):
    print(f"Plotting with tissue_type_name: {tissue_type_name}")  
    
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.scatter(
        all_coords[:, 0], 
        all_coords[:, 1], 
        all_coords[:, 2], 
        c='#C4D7FF', 
        s=0.8, 
        alpha=alpha_all, 
        label="All cells"
    )
    
    ax.scatter(
        highlight_coords[:, 0], 
        highlight_coords[:, 1], 
        highlight_coords[:, 2], 
        c='#1D1616', 
        s=0.9, 
        alpha=alpha_highlight, 
        label=tissue_type_name
    )
    
    ax.scatter(
        sender_receiver_data['Sender_x'], 
        sender_receiver_data['Sender_y'], 
        sender_receiver_data['Sender_z'], 
        c='#3D3BF3', 
        s=1.0,
        alpha=1.0, 
        label="Sender Cells"
    )
    ax.scatter(
        sender_receiver_data['Receiver_x'], 
        sender_receiver_data['Receiver_y'], 
        sender_receiver_data['Receiver_z'], 
        c='#FF2929', 
        s=1.0,
        alpha=1.0, 
        label="Receiver Cells"
    )
    
    for _, row in sender_receiver_data.iterrows():
        ax.plot(
            [row['Sender_x'], row['Receiver_x']],
            [row['Sender_y'], row['Receiver_y']],
            [row['Sender_z'], row['Receiver_z']],
            c='#B7B7B7', linestyle='-', linewidth=0.5, alpha=0.5
        )

    ax.set_xlabel("Spatial X Coordinates")
    ax.set_ylabel("Spatial Y Coordinates")
    ax.set_zlabel("Spatial Z Coordinates")
    ax.legend()
    
    ax.set_title(f"Pathway: {pathway_name}")
    
    if static_file:
        plt.savefig(static_file, bbox_inches='tight')
    plt.close()

def process_file_static(file_name, work_directory, args):
    print(f"Processing file: {file_name} with tissue type: {args.tissue_type}")  

    tissue_files_dir = os.path.join(work_directory, file_name, "tissue")
    pathway_dir = os.path.join(work_directory, file_name, "tissue/pathways")
    output_dir_static = os.path.join(work_directory, file_name, "plot/tissue/Plot_ligand_receptor_with_arrow/static")
    registered_csv = os.path.join(work_directory, file_name, f"{file_name}_registered.csv")
    
    registered_data = pd.read_csv(registered_csv, index_col=0)
    
    all_coords = registered_data[['spatial_registered_x', 'spatial_registered_y', 'spatial_registered_z']].dropna().to_numpy()
   
    os.makedirs(output_dir_static, exist_ok=True)
    
    for pathway_file_name in os.listdir(pathway_dir):
        if pathway_file_name.endswith('.csv') and args.tissue_type.lower() in pathway_file_name.lower():
            pathway_file_path = os.path.join(pathway_dir, pathway_file_name)
            
            sender_receiver_data = pd.read_csv(pathway_file_path)
            
            tissue_file_name = f"{args.tissue_type.strip()}.csv"
            tissue_file_path = os.path.join(tissue_files_dir, tissue_file_name)

            if os.path.exists(tissue_file_path):
                print(f"Found tissue file: {tissue_file_name}, loading tissue-specific cells...")
                
                tissue_data = pd.read_csv(tissue_file_path)
                print(f"Tissue data sample for {args.tissue_type}:")
                print(tissue_data.head())  
                
                matched_cells = registered_data.loc[registered_data.index.intersection(tissue_data['Unnamed: 0'])]

                if matched_cells.empty:
                    print(f"No matching cells found for tissue: {args.tissue_type}")
                else:
                    print(f"Number of matched cells: {matched_cells.shape[0]}")

                highlight_coords = matched_cells[['spatial_registered_x', 'spatial_registered_y', 'spatial_registered_z']].dropna().to_numpy()

                pathway_name = os.path.splitext(pathway_file_name)[0]
                static_file = os.path.join(output_dir_static, f"{pathway_name}.pdf")

                plot_cells_with_connections_static(
                    all_coords=all_coords,
                    highlight_coords=highlight_coords,
                    sender_receiver_data=sender_receiver_data,
                    tissue_type_name=args.tissue_type.strip(),
                    pathway_name=pathway_name,
                    static_file=static_file,
                    width=args.width,
                    height=args.height,
                    alpha_all=args.alpha_all,
                    alpha_highlight=args.alpha_highlight
                )
                print(f"Static plot saved for pathway {pathway_name} in {static_file}.")
            else:
                print(f"Tissue file {tissue_file_name} not found in {tissue_files_dir}")

def main_static(args):
    with open(args.filename_txt, 'r') as f:
        file_names = f.read().splitlines()

    for file_name in file_names:
        process_file_static(file_name, args.work_directory, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate static pathway plots for a specific tissue type in 3D.")
    parser.add_argument('--filename_txt', type=str, required=True, help="Path to the file containing the list of filenames.")
    parser.add_argument('--work_directory', type=str, required=True, help="Base working directory containing subfolders.")
    parser.add_argument('--tissue_type', type=str, required=True, help="Name of the tissue type to process (e.g., Brain, Heart).")
    parser.add_argument('--width', type=float, default=10, help="Width of the figure in inches.")
    parser.add_argument('--height', type=float, default=10, help="Height of the figure in inches.")
    parser.add_argument('--alpha_all', type=float, default=0.5, help="Transparency of all cells (0 to 1).")
    parser.add_argument('--alpha_highlight', type=float, default=1.0, help="Transparency of highlighted cells (0 to 1).")
    
    args = parser.parse_args()
    main_static(args)

