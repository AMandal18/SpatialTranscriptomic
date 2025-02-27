import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
from mpl_toolkits.mplot3d import Axes3D

def load_csv(file_path):
    """
    Load a CSV file into a pandas DataFrame.
    """
    return pd.read_csv(file_path)

def plot_3d_coordinates(data, filename=None, width=10, height=10, point_size=10, alpha=0.6):
    x = data['spatial_registered_x']
    y = data['spatial_registered_y']
    z = data['spatial_registered_z']
    
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='#3D3BF3', s=point_size, alpha=alpha)
    
    ax.set_xlabel('Spatial X Coordinates')
    ax.set_ylabel('Spatial Y Coordinates')
    ax.set_zlabel('Spatial Z Coordinates')
    
    if filename:
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        plt.savefig(filename, format='pdf')
        print(f"3D plot saved to {filename}")

    return ax

def main(args): 
    input_file = os.path.join(args.folderpath, f"{args.filename}_registered.csv")
    output_file = os.path.join(args.folderpath, "plot", "all_cells_plot_3D.pdf")
    
    if not os.path.exists(input_file):
        print(f"Input file not found: {input_file}")
        exit()
    
    data = load_csv(input_file)
    plot_3d_coordinates(
        data,
        filename=output_file,
        width=args.width,
        height=args.height,
        point_size=args.point_size,
        alpha=args.alpha
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot 3D spatial coordinates from a single CSV file.')
    parser.add_argument('--folderpath', type=str, required=True, help='Base folder path containing processed data.')
    parser.add_argument('--filename', type=str, required=True, help='Name of the subfolder containing the data (e.g., E8_5_rep1__Puck_201104_07).')
    parser.add_argument('--width', type=float, default=10, help='Width of the plot in inches.')
    parser.add_argument('--height', type=float, default=10, help='Height of the plot in inches.')
    parser.add_argument('--point_size', type=float, default=0.1, help='Size of the points in the scatter plot.')
    parser.add_argument('--alpha', type=float, default=0.6, help='Transparency of the points in the scatter plot.')
    args = parser.parse_args()
    main(args)

