import argparse
import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go

def plot_cells_with_connections_dynamic(
    all_coords,
    highlight_coords,
    sender_receiver_data,
    tissue_type_name,
    pathway_name,
    dynamic_file=None,
    alpha_all=0.8,
    alpha_highlight=0.9,
):
    print(f"Plotting with tissue_type_name: {tissue_type_name}")  
    
    fig_dynamic = go.Figure()
    
    fig_dynamic.add_trace(go.Scatter3d(
        x=all_coords[:, 0],
        y=all_coords[:, 1],
        z=all_coords[:, 2],
        mode='markers',
        marker=dict(size=0.9, color='#C4D7FF', opacity=alpha_all),
        name="All cells"
    ))
    
    fig_dynamic.add_trace(go.Scatter3d(
        x=highlight_coords[:, 0],
        y=highlight_coords[:, 1],
        z=highlight_coords[:, 2],
        mode='markers',
        marker=dict(size=0.9, color='#1D1616', opacity=alpha_highlight),
        name=tissue_type_name
    ))
    
    fig_dynamic.add_trace(go.Scatter3d(
        x=sender_receiver_data['Sender_x'],
        y=sender_receiver_data['Sender_y'],
        z=sender_receiver_data['Sender_z'],
        mode='markers',
        marker=dict(size=1, color='#3D3BF3', opacity=1.0),
        name="Sender Cells"
    ))
    
    fig_dynamic.add_trace(go.Scatter3d(
        x=sender_receiver_data['Receiver_x'],
        y=sender_receiver_data['Receiver_y'],
        z=sender_receiver_data['Receiver_z'],
        mode='markers',
        marker=dict(size=1, color='#FF2929', opacity=1.0),
        name="Receiver Cells"
    ))
    
    for _, row in sender_receiver_data.iterrows():
        fig_dynamic.add_trace(go.Scatter3d(
            x=[row['Sender_x'], row['Receiver_x']],
            y=[row['Sender_y'], row['Receiver_y']],
            z=[row['Sender_z'], row['Receiver_z']],
            mode='lines',
            line=dict(color='#B7B7B7', width=0.8),
            opacity=0.5,
            showlegend=False
        ))
    
    fig_dynamic.update_layout(
        title=f"Pathway: {pathway_name}",
        scene=dict(
            xaxis_title="Spatial X Coordinates",
            yaxis_title="Spatial Y Coordinates",
            zaxis_title="Spatial Z Coordinates"
        )
    )
    
    if dynamic_file:
        fig_dynamic.write_html(dynamic_file)

def process_file_dynamic(file_name, work_directory, args):
    print(f"Processing file: {file_name} with tissue type: {args.tissue_type}")  

    tissue_files_dir = os.path.join(work_directory, file_name, "tissue")
    pathway_dir = os.path.join(work_directory, file_name, "tissue/pathways")
    output_dir_dynamic = os.path.join(work_directory, file_name, "plot/tissue/Plot_ligand_receptor_with_arrow/dynamic")
    registered_csv = os.path.join(work_directory, file_name, f"{file_name}_registered.csv")
    
    registered_data = pd.read_csv(registered_csv, index_col=0)
    
    all_coords = registered_data[['spatial_registered_x', 'spatial_registered_y', 'spatial_registered_z']].dropna().to_numpy()
    
    os.makedirs(output_dir_dynamic, exist_ok=True)
    
    if not os.path.exists(pathway_dir):
        print(f"Pathway directory {pathway_dir} not found. Skipping {file_name}.")
        return
    
    for pathway_file_name in os.listdir(pathway_dir):
        if pathway_file_name.endswith('.csv') and args.tissue_type.lower() in pathway_file_name.lower():
            pathway_file_path = os.path.join(pathway_dir, pathway_file_name)
            
            sender_receiver_data = pd.read_csv(pathway_file_path)
            
            tissue_file_name = f"{args.tissue_type.strip()}.csv"
            tissue_file_path = os.path.join(tissue_files_dir, tissue_file_name)

            if os.path.exists(tissue_file_path):
                print(f"Found tissue file: {tissue_file_name}, loading tissue-specific cells...")
               
                tissue_data = pd.read_csv(tissue_file_path)
                
                matched_cells = registered_data.loc[registered_data.index.intersection(tissue_data['Unnamed: 0'])]

                if matched_cells.empty:
                    print(f"No matching cells found for tissue: {args.tissue_type}")
                else:
                    print(f"Number of matched cells: {matched_cells.shape[0]}")

                highlight_coords = matched_cells[['spatial_registered_x', 'spatial_registered_y', 'spatial_registered_z']].dropna().to_numpy()
                
                pathway_name = os.path.splitext(pathway_file_name)[0]
                dynamic_file = os.path.join(output_dir_dynamic, f"{pathway_name}.html")
                
                plot_cells_with_connections_dynamic(
                    all_coords=all_coords,
                    highlight_coords=highlight_coords,
                    sender_receiver_data=sender_receiver_data,
                    tissue_type_name=args.tissue_type.strip(),
                    pathway_name=pathway_name,
                    dynamic_file=dynamic_file,
                    alpha_all=args.alpha_all,
                    alpha_highlight=args.alpha_highlight
                )
                print(f"Dynamic plot saved for pathway {pathway_name} in {dynamic_file}.")
            else:
                print(f"Tissue file {tissue_file_name} not found in {tissue_files_dir}")

def main_dynamic(args):
    print(f"Reading filenames from: {args.filename_txt}")  
    
    with open(args.filename_txt, 'r') as f:
        file_names = f.read().splitlines()

    for file_name in file_names:
        process_file_dynamic(file_name, args.work_directory, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate dynamic pathway plots for a specific tissue type in 3D.")
    parser.add_argument('--filename_txt', type=str, required=True, help="Path to the file containing the list of filenames.")
    parser.add_argument('--work_directory', type=str, required=True, help="Base working directory containing subfolders.")
    parser.add_argument('--tissue_type', type=str, required=True, help="Name of the tissue type to process (e.g., Brain, Heart).")
    parser.add_argument('--alpha_all', type=float, default=0.5, help="Transparency of all cells (0 to 1).")
    parser.add_argument('--alpha_highlight', type=float, default=1.0, help="Transparency of highlighted cells (0 to 1).")

    args = parser.parse_args()
    main_dynamic(args)

