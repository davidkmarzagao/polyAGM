from polymer_regression_analysis import *
from smiles_to_grakel_graphs_chiral import *
# from weisfeiler_lehman_multiple_bonds import *
from weisfeiler_lehman_multiple_bond_search_instance import *
import grakel 
import pickle
from grakel.kernels import WeisfeilerLehman, VertexHistogram, WeisfeilerLehmanOptimalAssignment
import copy
import lime
import lime.lime_tabular
import matplotlib.pyplot as plt
from sklearn import metrics

import pyGPs
from pyGPs.Validation import valid
from pyGPs.GraphExtensions import graphUtil,graphKernels


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, ExpSineSquared
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold

from sklearn.feature_selection import RFE
import rdkit


from sklearn.model_selection import train_test_split
# import sklearn
from sklearn.model_selection import KFold
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import sklearn
import sklearn.ensemble 


from grakel.kernels import (
    GraphletSampling,
    RandomWalk,
    RandomWalkLabeled,
    ShortestPath,
    ShortestPathAttr,
    WeisfeilerLehman,
    WeisfeilerLehmanOptimalAssignment,
    NeighborhoodHash,
    PyramidMatch,
    SubgraphMatching,
    NeighborhoodSubgraphPairwiseDistance,
    LovaszTheta,
    SvmTheta,
    OddSth,
    Propagation,
    PropagationAttr,
    HadamardCode,
    MultiscaleLaplacian,
    VertexHistogram,
    EdgeHistogram,
    GraphHopper,
)

# N_JOBS is for paralelisation of kernel evaluation
N_JOBS = None

NORMALIZING_GRAPH_KERNELS = True


GRAKEL_KERNELS = {
    "SPath": lambda: ShortestPath(normalize=NORMALIZING_GRAPH_KERNELS),
    "VHist": lambda: VertexHistogram(normalize=NORMALIZING_GRAPH_KERNELS),
    "GSamp": lambda: GraphletSampling(normalize=NORMALIZING_GRAPH_KERNELS),
    "WL-1": lambda: WeisfeilerLehman(
        n_iter=1, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-2": lambda: WeisfeilerLehman(
        n_iter=2, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-3": lambda: WeisfeilerLehman(
        n_iter=3, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-4": lambda: WeisfeilerLehman(
        n_iter=4, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-5": lambda: WeisfeilerLehman(
        n_iter=5, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-6": lambda: WeisfeilerLehman(
        n_iter=6, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-1": lambda: WeisfeilerLehmanWithMultipleBonds(
        n_iter=1, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-2": lambda: WeisfeilerLehmanWithMultipleBonds(
        n_iter=2, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-3": lambda: WeisfeilerLehmanWithMultipleBonds(
        n_iter=3, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-4": lambda: WeisfeilerLehmanWithMultipleBonds(
        n_iter=4, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-5": lambda: WeisfeilerLehmanWithMultipleBonds(
        n_iter=5, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-6": lambda: WeisfeilerLehmanWithMultipleBonds(
    n_iter=6, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-7": lambda: WeisfeilerLehmanWithMultipleBonds(
    n_iter=7, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-10": lambda: WeisfeilerLehmanWithMultipleBonds(
    n_iter=10, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-20": lambda: WeisfeilerLehmanWithMultipleBonds(
    n_iter=20, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WLB-40": lambda: WeisfeilerLehmanWithMultipleBonds(
    n_iter=40, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "NH": lambda: NeighborhoodHash(
        n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "HC-1": lambda: HadamardCode(
        n_iter=1, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "HC-2": lambda: HadamardCode(
        n_iter=2, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "HC-3": lambda: HadamardCode(
        n_iter=3, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "HC-4": lambda: HadamardCode(
        n_iter=4, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "HC-5": lambda: HadamardCode(
        n_iter=5, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "HC-6": lambda: HadamardCode(
    n_iter=6, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    # There are no edge-labels
#     "NSPD": lambda: NeighborhoodSubgraphPairwiseDistance(
#         normalize=NORMALIZING_GRAPH_KERNEL
#     ),
    "WL-OA-1": lambda: WeisfeilerLehmanOptimalAssignment(
        n_iter=1, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-2": lambda: WeisfeilerLehmanOptimalAssignment(
        n_iter=2, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-3": lambda: WeisfeilerLehmanOptimalAssignment(
        n_iter=3, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-4": lambda: WeisfeilerLehmanOptimalAssignment(
        n_iter=4, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-5": lambda: WeisfeilerLehmanOptimalAssignment(
        n_iter=5, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-6": lambda: WeisfeilerLehmanOptimalAssignment(
    n_iter=6, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-7": lambda: WeisfeilerLehmanOptimalAssignment(
    n_iter=7, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-8": lambda: WeisfeilerLehmanOptimalAssignment(
    n_iter=8, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-20": lambda: WeisfeilerLehmanOptimalAssignment(
    n_iter=20, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "WL-OA-40": lambda: WeisfeilerLehmanOptimalAssignment(
    n_iter=40, n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "MS_Laplacian": lambda: MultiscaleLaplacian(
    n_jobs=N_JOBS, normalize=NORMALIZING_GRAPH_KERNELS
    ),
    "OddSth": lambda: OddSth(normalize=NORMALIZING_GRAPH_KERNELS),
#     "RUnits": lambda: create_repeat_units_feature_vectors(),
    "RUnits": lambda: (),
}



from sklearn.decomposition import PCA
from sklearn import linear_model
import pysmiles
import csv  

def restrict_list_of_graphs_based_on_non_null_targets(rows_of_block_ratios, graphs, target, rows_with_null_block_ratios):
    
    provisory_y = []
    final_y = []
    final_G = [] # Where the graphs associated to non-null entries y will be added to. 
    #  This checks y entries for the non-null block_ration entries. The y entries can be null for now. 
    for counter in range(len(rows_with_null_block_ratios)):
        if rows_with_null_block_ratios[counter]:
            continue
        provisory_y.append(float((target[counter])))
    #  This now checks entries of y against the non-null entries of block ratio. 
    #  A new graphs vector will be created based on the non-null values of y
    y_is_null = np.isnan(provisory_y)
    for i in range(len(y_is_null)):
        if y_is_null[i]:
            continue
        final_y.append(provisory_y[i])
        final_G.append(graphs[i])
        
    return final_G, final_y


def create_repeat_units_feature_vectors():
    smiles = df['Big_Smile']
    block_ratios = df['Block_ratio']
    repeat_unit_dict = generate_dictionary_of_repeat_units(smiles)
    return from_smiles_to_repeat_units(smiles, block_ratios, rows_with_null_block_ratios, repeat_unit_dict)

def generate_explanations_and_subgraphs(exp, gk, importance_order, i, d, index_of_tests, index_of_train):
    if importance_order == 0:  
        exp.show_in_notebook(show_table=True)
        print(exp.as_list())
    
    print(exp.as_list()[importance_order])
    pattern, change = split_explanation(exp.as_list()[importance_order])
    print(pattern)
    print(change)
    found_in_test_graphs = False
    found_in_train_graphs = False
    graph_where_pattern_was_found = None
    print("index of tests: ", index_of_tests)
    for test_graph in index_of_tests:
        if found_in_test_graphs:
            break
        for node in gk._patterns[d][test_graph]:
#             print(pattern)
#             print(gk._patterns[d][test_graph][node])
            if pattern == gk._patterns[d][test_graph][node]:
                print(pattern)
                print(gk._patterns[d][test_graph][node])
                print('Found')
                found_in_test_graphs = True
                graph_where_pattern_was_found = test_graph
                print(node)
                print(test_graph)
                break
    if not found_in_test_graphs:
        print("pattern not found in test graphs")
    for train_graph in index_of_train:
        if found_in_test_graphs:
            break
        if found_in_train_graphs:
            break
        for node in gk._patterns[d][train_graph]:
#             print(pattern)
#             print(gk._patterns[d][test_graph][node])
            if pattern == gk._patterns[d][train_graph][node]:
                print(pattern)
                print(gk._patterns[d][train_graph][node])
#                 print('Found')
                found_in_train_graphs = True
                graph_where_pattern_was_found = train_graph
                print(node)
                print(train_graph)
                break
    if not found_in_train_graphs:
        print("pattern not found in train graphs")
    return node, graph_where_pattern_was_found, pattern, found_in_test_graphs, change

def generate_networkx_subgraph_from_grakel_graph_node_and_depth(gr, depth, starting_node):
    subgraph = nx.Graph()
    # depth = 2
    queue = []
    # queue.append(13)
    queue.append(starting_node)
    # print(gr[0])
    for d in range(depth):
        n_nodes_to_add = len(queue) # To make sure we don't add neighbours of neighbours at this depth. 
        for i in range(n_nodes_to_add):
            new_node = queue.pop(0)
    #         print(gr[0][queue.pop(0)])
            if new_node not in list(subgraph.nodes):
                subgraph.add_node(new_node, element= gr[1][new_node])
#             print(gr[0][new_node])
            for neigh, weight in gr[0][new_node].items():
                if neigh not in list(subgraph.nodes):
                    subgraph.add_node(neigh, element= gr[1][neigh])
                    queue.append(neigh)
                subgraph.add_edge(new_node, neigh, order = weight)
    return subgraph

def defining_learning_method(targ, l_met, weight_y_data_train, weight_y_data_train_tg_tg2, y_var_data_train):
    if l_met == 'rf':
        met = sklearn.ensemble.RandomForestRegressor(n_estimators=1000)
    if l_met == 'ard':
        met = linear_model.ARDRegression()
    if l_met == 'bayes':
        met = linear_model.BayesianRidge()
    if l_met == 'gp-scikit':
#         kernel = 5*DotProduct(1, sigma_0_bounds = (10  **(-2), 10 **1 ))+ WhiteKernel(0.5, noise_level_bounds=(10  **(-2), 10 **(-1))) +  RBF(length_scale_bounds = (10  **(-2), 3*10 **1 )) #+ Matern()
#         kernel = 5*DotProduct(1, sigma_0_bounds = (10  **(-3), 10 **3 ))+ WhiteKernel(0.5, noise_level_bounds=(10  **(-4), 10 **(-1))) +  RBF(length_scale_bounds = (10  **(-3), 3*10 **1 )) #+ Matern()
        kernel = DotProduct(1, sigma_0_bounds = (10  **(-2), 10 **3 ))+ WhiteKernel(0.5, noise_level_bounds=(10  **(-4), 10 **(-1))) +  RBF(length_scale_bounds = (10  **(-3), 3*10 **1 )) #+ Matern()
#         kernel = 1*DotProduct(1, sigma_0_bounds = (10  **(-3), 10 **3 )) +  RBF(length_scale_bounds = (10  **(-3), 3*10 **1 )) # + Matern()
        # Fixing noise for Tg and Tg2 to a fixed value:
        if t == 'Tg' or t == 'Tg2':
            y_var_data_train = np.where(y_var_data_train == 0, weight_y_data_train_tg_tg2, y_var_data_train)
        elif t in ['E_(MPa)', 'σbreak (MPa)', 'εbreak (pct)']:
            y_var_data_train = weight_y_data_train * y_var_data_train
        met = GaussianProcessRegressor(kernel=kernel, random_state=0, alpha = np.square(y_var_data_train), normalize_y=True, n_restarts_optimizer=5)
    return met
    
def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

def add_row_to_csv_file_of_smiles(new_row, file_name, is_first_row):
    with open(file_name, 'a') as f:
        # check if file em empty, if it is, add title columns
        writer = csv.writer(f)
        if is_first_row: 
            writer.writerow(['change', 'smiles'])
        writer.writerow(new_row)
    




def plot_predictions_against_real_values_with_hovering(x, y, error_variance, experimental_error_variance, ker, t, index_of_test_data, pol_names, pol_numbers, pol_colours, error_bars = False, x_error_bars = True):
    # x is the real value. y is the predicted 
    
    names = [pol_names[int] for int in index_of_test_data]
    numbers = [str(pol_numbers[int]) for int in index_of_test_data]
    colours = [pol_colours[int] for int in index_of_test_data] 

    fig = plt.figure()
#     ax = fig.add_axes([0.1, 0.1, 0.3, 0.7])
#     ax = fig.add_axes([0.1, 0.1, 0.3, 0.5])
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#     colours_to_use = []
#     for i in range(len(c)):
#         colours_to_use.append(colours[colour_names[c[i]]])
#     fig,ax = plt.subplots()
    
#     ax = fig.add_axes([0.1, 0.1, 0.5, 0.75])
#     fig,ax = px.subplots()
    for i in range(len(x)):
#         sc = ax.scatter(x[i], y[i], c=colours[colour_names[c[i]]], alpha=0.5, 
        sc = ax.scatter(x[i], y[i], c=colours[i], alpha=0.5, 
#                          cmap=cmap,  
                         s=80, 
                         label=names[i]
                        )
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1, 0.5), ncol = 1)
    
    
    sc_2 = ax.scatter(x, y, c=colours, alpha=0.5, 
                     s=80, 
                    )

    int(max(x)) + 1
    print(mean_squared_error(x, y))
    if t == 'σbreak (MPa)':
        x_max = int(max(x)) + 3
        plt.plot(range(int(min(x) - 1), int(max(x)) + 4), range(int(min(x) - 1), int(max(x)) + 4), color="black")
    else:
        x_max = int(max(x))
    plt.plot(range(int(min(x) - 1), x_max + 1), range(int(min(x) - 1), x_max + 1), color="black")
    
    plt.text(int(min(x) + 1), x_max, "RMSE: {} \nR$^2$: {}".format(round(mean_squared_error(x, y, squared=False), 2), round(r2_score(x, y), 2)), verticalalignment='top', fontsize=10)

        
        
    plt.xlabel("Experimental " + pretty_target_names[t], fontsize=10)
    plt.ylabel("Predicted " + pretty_target_names[t], fontsize=10)
#     plt.title("Experimental vs Predicted " + pretty_target_names[t] + " with " + ker)
    add_title = False
    if add_title:
        plt.title("Experimental vs Predicted " + pretty_target_names[t])
    


    
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)
    

    
    def update_annot(ind):

        pos = sc_2.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}, {}".format(" ".join([numbers[n] for n in ind["ind"]]), 
                                " ".join([names[n] for n in ind["ind"]]))
        annot.set_text(text)
#         annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
#         annot.get_bbox_patch().set_facecolor(cmap(c[ind["ind"][0]]))
        annot.get_bbox_patch().set_facecolor(colours[ind["ind"][0]])
        annot.get_bbox_patch().set_alpha(0.4)
    
   

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc_2.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()
    
    
    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.axis('square')
    
    plt.show()
    
#     print(error_variance)
    if error_bars:
        d= np.sqrt(error_variance)
        #     print(d)
        plt.errorbar(x,y,yerr=d, linestyle="None")
        
    # This is for the experimental error reported in the papers. 
    if x_error_bars:
#         error_x = np.sqrt(experimental_error_variance)
        error_x = experimental_error_variance
#         plt.errorbar(x,y, xerr=error_x, ecolor='gray', linestyle="None")
        markers, caps, bars = plt.errorbar(x,y, xerr=error_x,  color="b", linestyle="None")
        [bar.set_alpha(0.5) for bar in bars]
#     print("Writing figure in HTML")

    fig.savefig("plot_in_svg.svg", format="svg", dpi=2400)
    plt.show()
#     plt.tight_layout()


def split_explanation(exp):
    average_change = exp[1]
    strng = exp[0]
    fingerprint = ''
    one_sided_inequality = True
    if exp[0].count('<') == 2:
        one_sided_inequality = False
        ineq = '<'
    else:
        if exp[0].count('>') == 2:
            one_sided_inequality = False
            ineq = '>'
        else:
            one_sided_inequality = True
    if not one_sided_inequality:
        fingerprint = strng[(strng.find(ineq) + 2):]
        print(fingerprint)
        fingerprint = fingerprint[:(fingerprint.find(ineq) - 1)]
        print(fingerprint)
#         exp[0] = exp[0][exp[0].find(ineq):]
#         [s.find('egg'):]
#         print(exp[0])
#         exp[0] = exp[0][2:]
#         print(exp[0])
    else:
        for char in exp[0]:
#         print(char)
            if char not in ['<', '>', '=']:
                fingerprint = fingerprint + str(char)
            else: 
                fingerprint = fingerprint[:-1]  # There was an extra space char at the end that we want to delete
                break
    return fingerprint, average_change



pretty_target_names = {
        'Tg' : 'Lower $T_g$',
         'Tg2': 'Upper $T_g$',
         'σbreak (MPa)' : '$σ_{break}$ (MPa)',
         'εbreak (pct)' : '$ε_{break}$ (%)',
}




