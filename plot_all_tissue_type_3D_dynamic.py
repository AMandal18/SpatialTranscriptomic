import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
import re

def load_tissue_colors(color_csv):
    """
    Loads tissue colors from a CSV file and ensures no whitespace issues.
    """
    tissue_color_df = pd.read_csv(color_csv, header=None, names=["tissue_type", "color"])
    tissue_colors = {k.strip(): v.strip() for k, v in zip(tissue_color_df["tissue_type"], tissue_color_df["color"])}
    return tissue_colors


def standardize_tissue_name(tissue_name):
    """
    Standardizes tissue names by replacing underscores with spaces and removing special characters.
    """
    tissue_name = tissue_name.replace("_", " ")  
    tissue_name = re.sub(r"\(.*?\)", "", tissue_name).strip()  
    return tissue_name


def plot_all_tissues_3d(all_tissue_coords, tissue_colors, title="All Tissue Types in 3D Space"):
    """
    Creates a 3D scatter plot with different colors for each tissue type.
    """
    fig = go.Figure()

    for tissue, coords in all_tissue_coords.items():
        if coords.shape[0] > 0:
            standardized_tissue = standardize_tissue_name(tissue.replace(".csv", ""))  
            color = tissue_colors.get(standardized_tissue, "gray")

            print(f"Plotting {standardized_tissue} with color {color}")  

            fig.add_trace(
                go.Scatter3d(
                    x=coords[:, 0],
                    y=coords[:, 1],
                    z=coords[:, 2],
                    mode='markers',
                    marker=dict(size=3, color=color, opacity=0.8),
                    name=standardized_tissue
                )
            )

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='Spatial X Coordinates',
            yaxis_title='Spatial Y Coordinates',
            zaxis_title='Spatial Z Coordinates'
        ),
        template='plotly_white'
    )
    return fig


def process_single_file(file_name, work_directory, tissue_colors, tissue_list_path):
    tissue_files_dir = os.path.join(work_directory, file_name, "tissue")
    registered_csv = os.path.join(work_directory, file_name, f"{file_name}_registered.csv")

    if not os.path.exists(registered_csv):
        print(f"Warning: {registered_csv} not found. Skipping {file_name}.")
        return

    registered_data = pd.read_csv(registered_csv, index_col=0)

    with open(tissue_list_path, 'r') as f:
        tissue_files = [line.strip().replace(' ', '_') + '.csv' for line in f]

    all_tissue_coords = {tissue: [] for tissue in tissue_files}

    for tissue_file in tissue_files:
        tissue_file_path = os.path.join(tissue_files_dir, tissue_file)
        if os.path.exists(tissue_file_path):
            tissue_data = pd.read_csv(tissue_file_path)
            matched_cells = registered_data.loc[registered_data.index.intersection(tissue_data['Unnamed: 0'])]
            highlight_coords = matched_cells[['spatial_registered_x', 'spatial_registered_y', 'spatial_registered_z']].dropna().to_numpy()

            if highlight_coords.shape[0] > 0:
                all_tissue_coords[tissue_file].append(highlight_coords)

    filtered_tissue_coords = {
        tissue: np.vstack(coords) for tissue, coords in all_tissue_coords.items() if len(coords) > 0
    }

    if len(filtered_tissue_coords) == 0:
        print(f"No tissue data found for {file_name}, skipping plot generation.")
        return
    
    combined_fig = plot_all_tissues_3d(filtered_tissue_coords, tissue_colors)
    
    output_dir = os.path.join(work_directory, file_name, "plot", "tissue")
    os.makedirs(output_dir, exist_ok=True)  
    
    output_html = os.path.join(output_dir, "all_tissues_3D_plot.html")
    combined_fig.write_html(output_html)
    print(f"Saved: {output_html}")


def main(args):
    tissue_colors = load_tissue_colors(args.color_csv)
    
    with open(args.filename_txt, 'r') as f:
        file_names = f.read().splitlines()
    
    for file_name in file_names:
        print(f"\nProcessing {file_name}...")
        process_single_file(file_name, args.work_directory, tissue_colors, args.tissue_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot all tissue types in separate 3D plots.")
    parser.add_argument('--filename_txt', type=str, required=True, help="Path to the file containing the list of folder names.")
    parser.add_argument('--work_directory', type=str, required=True, help="Base working directory containing subfolders.")
    parser.add_argument('--color_csv', type=str, required=True, help="CSV file containing tissue types and colors.")
    parser.add_argument('--tissue_list', type=str, required=True, help="File containing list of tissue names.")
    args = parser.parse_args()
    main(args)

