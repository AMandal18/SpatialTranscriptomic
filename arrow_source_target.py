import argparse
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def plot_arrow_network(csv_file, source_cell, target_cell, output_file):
    df = pd.read_csv(csv_file)
    
    required_columns = {'source', 'ligand', 'receptor', 'target', 'score'}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"CSV file must contain the following columns: {required_columns}")
    
    filtered_df = df[(df['source'] == source_cell) & (df['target'] == target_cell)]

    if filtered_df.empty:
        print(f"No data found for source '{source_cell}' and target '{target_cell}'.")
        return
    
    G = nx.DiGraph()
    
    for _, row in filtered_df.iterrows():
        G.add_node(row['source'], layer=0, shape="rect")  
        G.add_node(row['ligand'], layer=1, shape="circle")  
        G.add_node(row['receptor'], layer=2, shape="circle")  
        G.add_node(row['target'], layer=3, shape="rect")  
        
        G.add_edge(row['source'], row['ligand'], weight=row['score'])
        G.add_edge(row['ligand'], row['receptor'], weight=row['score'])
        G.add_edge(row['receptor'], row['target'], weight=row['score'])
    
    pos = {}
    node_layers = {
        0: [source_cell], 
        1: list(filtered_df['ligand'].unique()), 
        2: list(filtered_df['receptor'].unique()), 
        3: [target_cell]
    }
    
    for layer, nodes in node_layers.items():
        for i, node in enumerate(nodes):
            if layer == 1:  
                pos[node] = (layer + 0.2, -i * 0.8)  
            elif layer == 2:  
                pos[node] = (layer - 0.2, -i * 0.8)  
            else:
                pos[node] = (layer, -i * 0.8)  
    
    ligands = node_layers[1]
    receptors = node_layers[2]

    middle_ligand_y = pos[ligands[len(ligands) // 2]][1] if ligands else 0
    middle_receptor_y = pos[receptors[len(receptors) // 2]][1] if receptors else 0
    
    pos[source_cell] = (0 + 1.05, middle_ligand_y)
    pos[target_cell] = (3 - 1.05, middle_receptor_y)
    
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    
    min_width, max_width = 0.5, 5 
    if max(edge_weights) == min(edge_weights):  
        norm_weights = [max_width] * len(edge_weights)  
    else:
        norm_weights = [(max_width - min_width) * (w - min(edge_weights)) / (max(edge_weights) - min(edge_weights)) + min_width for w in edge_weights]
    
    plt.figure(figsize=(6, 5))
    
    nx.draw(G, pos, with_labels=False, node_color='white', edge_color='black', node_size=2000, font_size=10, font_weight='bold', arrows=True)
    nx.draw_networkx_edges(G, pos, width=norm_weights, edge_color='black', arrows=True)
    
    for node, (x, y) in pos.items():
        if node == source_cell:  
            formatted_label = '\n'.join(source_cell.split(' ')) if ' ' in source_cell else source_cell
            plt.text(x - 0.01, y, formatted_label, verticalalignment='center', horizontalalignment='center',
                     fontsize=15, fontweight='bold', family='serif', rotation=90, bbox=dict(facecolor="white", edgecolor="black"),
                     color='blue')
        elif node == target_cell:  
            formatted_label = '\n'.join(target_cell.split(' ')) if ' ' in target_cell else target_cell
            plt.text(x + 0.03, y, formatted_label, verticalalignment='center', horizontalalignment='center',
                     fontsize=15, fontweight='bold', family='serif', rotation=90, bbox=dict(facecolor="white", edgecolor="black"),
                     color='red')
        else:
            color = 'blue' if node in node_layers[1] else 'red'
            plt.text(x, y, node, fontsize=12, fontweight='bold', family='serif', ha='center', color=color)
    
    for u, v in G.edges():
        if u == source_cell:
            pos[v] = (pos[v][0], middle_ligand_y)  
        if v == target_cell:
            pos[u] = (pos[u][0], middle_receptor_y)  

    
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

    print(f"Arrow plot saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an arrow plot matching the reference style.")
    parser.add_argument("--csv_file", type=str, required=True, help="Input CSV file containing source, ligand, receptor, target, and score.")
    parser.add_argument("--source", type=str, required=True, help="Source cell type.")
    parser.add_argument("--target", type=str, required=True, help="Target cell type.")
    parser.add_argument("--output_file", type=str, required=True, help="Output filename (e.g., output.pdf or output.png).")

    args = parser.parse_args()

    plot_arrow_network(args.csv_file, args.source, args.target, args.output_file)

