import itertools
import numpy as np
import networkx  as nx
import pickle
import json
import matplotlib.pyplot as plt
from matplotlib import gridspec

import argparse
import warnings
#import argcomplete
 

warnings.filterwarnings("ignore") 

parser = argparse.ArgumentParser(description="Plot any pharmacophore json file in graph formats")

## input coordinate files
parser.add_argument('-j', dest='pharmacophore_file', type=str, default='pharma.json',help=' pharmacophore file')
 
parser.add_argument('-o', dest='out_filename', type=str, default='./test.png',help='output figure')

parser.add_argument('-set_label', dest='set_label', type=int, default=1,help='Whether to set labels and scale bar or not')
#############################################################
args                   = parser.parse_args()
pharmacophore_file     = args.pharmacophore_file
out_filename           = args.out_filename
set_label              = args.set_label

with open(pharmacophore_file, 'r') as file: pharmacophore_data = json.load(file)


# Node color dictionary
node_colors =dict(Aromatic="orange", HydrogenDonor="magenta", HydrogenAcceptor="green", Hydrophobic="cyan",
                  PositiveIon="red", ExclusionSphere="gray", NegativeIon="blue")

def get_distance (point1, point2):
    return np.linalg.norm(np.array(point1) - np.array(point2))

def get_colors_for_edges (graph, colormap='coolwarm'):
    cmap = plt.colormaps[colormap]
    Edge_weights = np.array([weight for i, j, weight in graph.edges(data="distance")])
    min_dist, max_dist = min (Edge_weights), max (Edge_weights)
    Edge_weights  /= max(Edge_weights) 
    edge_color=[cmap(weight) for weight in Edge_weights]

    return edge_color, min_dist, max_dist

def get_node_color (graph):
    colors = [x[1] for x in list(graph.nodes (data="color"))]
    return colors
########################################################


# create graph from pharmacophore

G = nx.Graph()
Exclusion_spheres = []
index = 0

for sites in pharmacophore_data ["points"]:
    
    if sites ['name'] != 'ExclusionSphere':

        index += 1
        # Value should be the importance of a node or something. Like energy 
        G.add_node (index,
                    name= sites ['name'] , 
                    position=np.array ([sites ['x'], sites ['y'], sites ['z']]), 
                    radius = sites ['radius'],
                    size = sites['size'],
                    requirement = sites['requirement'],
                    enabled=sites['enabled'],
                    svector=sites['svector'],
                    vector=sites['vector'],
                    minsize = sites ['minsize'],
                    maxsize=sites ['maxsize'],
                    vector_on= sites['vector_on'],
                    color=node_colors [sites ['name']], 
                    score=sites ['score'],
                    chainID=sites['chainID'])

    else: Exclusion_spheres.append (sites)
        
possible_combinations = list (itertools.combinations ( list(G.nodes), 2))


for i,j in possible_combinations:
    
    score=G.nodes[i]['score']['normed_score'] + G.nodes[j]['score']['normed_score'] # score should be something to denote the importance of an edge or node
    distance = get_distance (G.nodes[i]['position'], G.nodes[j]['position'])
    G.add_edge (i,j, distance= distance, score=score)

  # initialize graphs
 
 # defint the node colors
node_colors = get_node_color (G)
edge_color , min_dist, max_dist = get_colors_for_edges (G)
colormap='coolwarm'
cmap = plt.colormaps[colormap]
#####################################################################3
node_edge_colors = ['k' if G.nodes[n]['chainID'] == 'A' else 'red' for n in G.nodes()]  # Different edge colors


#### Plots
if set_label:
    fig    =     plt.figure(figsize=(5.5,4))
else:
    fig    =     plt.figure(figsize=(4,4))
    
gs     =     fig.add_gridspec(nrows=1, ncols=1,width_ratios=[1], wspace=0.2, hspace=0.5)
ax     =     fig.add_subplot(gs[0, 0])

nx.draw (G, pos=nx.circular_layout(G), node_color=node_colors, with_labels=True, ax=ax, edge_color=edge_color,
        linewidths=2,width=1.5,edgecolors=node_edge_colors)

# plot the colormap
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_dist, vmax=max_dist))
sm.set_array([])
if set_label:
    cbar = plt.colorbar(sm, ax=ax,aspect=18, shrink=0.85,location='left')
    cbar.set_label('Node distance [$\mathrm{\AA}$]')

#path = "/mnt/second/chen-zacharias-data/apo/gromacs-traj/6ygj/new-search-feb-2025/with-aromatic/"
#plt.savefig (path+"./parent_graph.png", dpi=300, bbox_inches="tight") 
 
site_type = set()
for _, node_attr in G.nodes (data=True):
    site = node_attr  ['name']
    site_type.add (site)

  # initialize graphs
# Node color dictionary
node_colors =dict(Aromatic="orange", HydrogenDonor="magenta", HydrogenAcceptor="green", Hydrophobic="cyan",
                  PositiveIon="red", ExclusionSphere="gray", NegativeIon="blue")

if set_label:
    pos = [(0.1,1.07), (0.7,1.07),
          (0.1,1.14), (0.7,1.14),  (0.1,1.21)]
    for i, s in enumerate (site_type):
        print (s)
        ax.annotate (s, xycoords="axes fraction", xy=pos[i], color=node_colors [s], fontsize=12)
#out_path = "/mnt/second/chen-zacharias-data/apo/gromacs-traj/5ad3/new-search-feb-2025/figures/"
 
plt.savefig (out_filename, bbox_inches="tight", dpi=600);
