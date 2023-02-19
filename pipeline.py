df = pd.read_excel('database.xlsx')

#  We clean the data but do not change the slashed in block ratio. We also create another df with variances (zero if not reported)
df, df_var = cleanDataWithVariances(df, ['±', '/'], ['Block_ratio', 'Hard_block', 'Soft_block', 'ignore_row', 'sec', 'sterochemistry' ] )


df = df.infer_objects()
pd.set_option('display.max_rows', None)
pd.set_option('mode.chained_assignment', None)

#  possible to also print df_var
#df_var
#  printing df
#df

#  Identifying rows without block ratios. This will later be used to match the appropriate target property with its related polymer
rows_with_null_block_ratios = df['Block_ratio'].isnull()
# print(rows_with_null_block_ratios)


rows_to_ignore = []

num_entries = len(df['Monomers'])
for row in range(num_entries):
    if df['ignore_row'][row] == 1:
        rows_to_ignore.append(row)

# print(rows_to_ignore)     

rows_to_also_ignore = [27,64, 65,66,67,68, 86, 89, 90]

rows_to_ignore = rows_to_ignore + rows_to_also_ignore

print(rows_to_ignore)


for r_to_ignore in rows_to_ignore:
    rows_with_null_block_ratios[r_to_ignore] = True
# print(rows_with_null_block_ratios)    

smiles = df['Big_Smile']

# smiles = [strig]
# block_ratios = ['1:1:1']

block_ratios = df['Block_ratio']
print(block_ratios)


rows_with_stereo = [(df['stereochemistry'][row] == 1) for row in range(num_entries)]
# rows_with_stereo

CHIRALITY = True
HYDROGENS = False

full_pols = from_smiles_to_networkx(smiles, block_ratios, rows_with_null_block_ratios, rows_with_stereo, use_chiral = CHIRALITY, explicit_H = HYDROGENS)



# Gettings names of polymers with no null block ratio
# names = df['Monomers']
names = df['Hard_block'].astype(str) +"/"+ df["Soft_block"]
print(names)
colours = ['black', 'lightcoral', 'sandybrown', 'indianred', 'brown', 'chocolate', 'gold', 'fuchsia', 'darkkhaki', 'olive', 'steelblue', 'darkorange','darkolivegreen','deepskyblue', 'tan', 'moccasin', 'orange','red', 'gold', 'darkkhaki', 'olive', 'forestgreen', 'skyblue', 'rebeccapurple', 'yellow', 'yellowgreen', 'greenyellow', 'chartreuse', 'lawngreen', 'honeydew', 'darkseagreen', 'palegreen', 'lightgreen', 'forestgreen', 'skyblue', 'rebeccapurple', 'limegreen', 'darkgreen', 'green', 'lime', 'seagreen', 'mediumseagreen', 'springgreen', 'mintcream', 'mediumspringgreen', 'mediumaquamarine', 'aquamarine', 'turquoise', 'lightseagreen', 'mediumturquoise', 'azure', 'lightcyan', ]
# colours = ['black', 'lightcoral', 'indianred', 'brown', 'firebrick', 'maroon', 'darkred', 'red', 'mistyrose', 'salmon', 'tomato', 'darksalmon', 'coral', 'orangered', 'lightsalmon', 'sienna', 'seashell', 'chocolate', 'saddlebrown', 'sandybrown', 'peachpuff', 'peru', 'linen', 'bisque', 'darkorange', 'burlywood', 'tan', 'papayawhip', 'moccasin', 'orange', 'wheat', 'oldlace', 'darkgoldenrod', 'goldenrod', 'cornsilk', 'gold', 'lemonchiffon', 'darkkhaki', 'lightyellow', 'lightgoldenrodyellow', 'olive', 'yellow', 'olivedrab', 'yellowgreen', 'darkolivegreen', 'greenyellow', 'chartreuse', 'lawngreen', 'honeydew', 'darkseagreen', 'palegreen', 'lightgreen', 'forestgreen', 'limegreen', 'darkgreen', 'green', 'lime', 'seagreen', 'mediumseagreen', 'springgreen', 'mintcream', 'mediumspringgreen', 'mediumaquamarine', 'aquamarine', 'turquoise', 'lightseagreen', 'mediumturquoise', 'azure', 'lightcyan', 'paleturquoise', 'darkslategray', 'darkslategrey', 'teal', 'darkcyan', 'aqua', 'cyan', 'darkturquoise', 'cadetblue', 'powderblue', 'lightblue', 'deepskyblue', 'skyblue', 'lightskyblue', 'steelblue', 'aliceblue', 'dodgerblue', 'lightslategray', 'lightslategrey', 'slategray', 'slategrey', 'lightsteelblue', 'cornflowerblue', 'royalblue', 'ghostwhite', 'lavender', 'midnightblue', 'navy', 'darkblue', 'mediumblue', 'blue', 'slateblue', 'darkslateblue', 'mediumslateblue', 'mediumpurple', 'rebeccapurple', 'blueviolet', 'indigo', 'darkorchid', 'darkviolet', 'mediumorchid', 'thistle', 'plum', 'violet', 'purple', 'darkmagenta', 'fuchsia', 'magenta', 'orchid', 'mediumvioletred', 'deeppink', 'hotpink', 'lavenderblush', 'palevioletred', 'crimson', 'pink', 'lightpink']


name_to_colour_dict = dict()
counter = 0
for name in names:
    if name not in name_to_colour_dict.keys():
        name_to_colour_dict[name] = colours[counter]
        counter += 1 
print(name_to_colour_dict)
    
polymer_names = []
polymer_numbers = []
polymer_colours = []
for counter in range(len(names)):
    if len(rows_with_null_block_ratios) != 0:
        if rows_with_null_block_ratios[counter]:
            continue
    polymer_names.append(names[counter]) 
    polymer_numbers.append(counter) 
    polymer_colours.append(name_to_colour_dict[names[counter]])

print(polymer_numbers)
# polymer_names

#  DF with structures of monomers and smiles. 
monomers_df = pd.read_excel('monomers with smiles and structures(1).xlsx')
monomers_df

# Creating dict of whether a monomer is part of the hard block
is_hard_block = dict(zip(monomers_df.Monomer, monomers_df.Hard_Block))
print(is_hard_block)

# Creating dict of polymer SMILEs
smiles_of_monomer = dict(zip(monomers_df.Monomer, monomers_df.SMILES))
print(smiles_of_monomer)

# Creating dict of polymer molar mass
molar_mass_of_monomer = dict(zip(monomers_df.Monomer, monomers_df.Molar_mass))
print(molar_mass_of_monomer)



from grakel.utils import graph_from_networkx

#  Transforming networkX graphs into grakel graphs. 
grakel_pols = graph_from_networkx(full_pols, node_labels_tag='element', edge_weight_tag ='order', as_Graph = False)


#  Getting a list of them rather than an itearable. 

# G = grakel_pols
Graphs = []
for g in grakel_pols:
    Graphs.append(g)
Graphs


def run_kernels_against_target(df, 
                               target, 
                               target_var, 
                               graphs, 
                               number_of_runs, 
                               kernels = ["WL-4"], 
                               learning_met = 'rf', 
                               error = 'mean_sq', 
                               normalise = True, 
                               explicit_features = False,
                               explanation = False,
                               graph_being_explained = 70,
                               name_file_explanations = 'tg2_importances.csv',
                               extra_tests = [],
                               visualisation = False,
                               show_error_bars = False,
                              training_error = False,
                               weight_y_var_train = 1,
                               weight_y_var_train_for_tg_and_tg2 = 1,
                               truncate_to_zero_if_neg = True,
                              log_y = True, 
                              printing_results = False,
                              with_shortest_paths = False,
                              n_shortest_paths_features = 50,
                              RFE_selection = False,
                              RFE_selection_n_features = 50,
                              covar_feature_selection = False,
                               concatenate_depths = False):
    means = []
    stds = []
    size_of_test = 0.1
    runs_for_plot = round((1 / size_of_test) + 1, 0)
    
    G, y = restrict_list_of_graphs_based_on_non_null_targets(rows_with_null_block_ratios, graphs, target, rows_with_null_block_ratios)
    G, y_var = restrict_list_of_graphs_based_on_non_null_targets(rows_with_null_block_ratios, graphs, target_var,rows_with_null_block_ratios)
    if len(y) != len(y_var):
        print("Lenght of y and y_var are not the same")
    if len(G) != len(y):
        print("Lenght of G and y are not the same")    
    # Same function also works for restricting the list of polymer names and numbers. 
    names_of_polymers, y = restrict_list_of_graphs_based_on_non_null_targets(rows_with_null_block_ratios, polymer_names, target, rows_with_null_block_ratios)
    numbers_of_polymers, y = restrict_list_of_graphs_based_on_non_null_targets(rows_with_null_block_ratios, polymer_numbers, target, rows_with_null_block_ratios)
    colours_of_polymers, y = restrict_list_of_graphs_based_on_non_null_targets(rows_with_null_block_ratios, polymer_colours, target, rows_with_null_block_ratios) 
    if with_shortest_paths:
        features_shortest_paths, y = restrict_list_of_graphs_based_on_non_null_targets(rows_with_null_block_ratios, features_sp[:, 0:n_shortest_paths_features], target, rows_with_null_block_ratios)
    if len(y) != len(y_var):
        print("Lenght of y and y_var are not the same")
    if len(G) != len(y):
        print("Lenght of G and y are not the same")    
    
    if len(extra_tests) != 0: 
        G_extra_tests = extra_tests
            
    print("> Number of data points for", t, "is", len(y))
    print("> With shortest paths = ", with_shortest_paths)
    print("> With concatenate depth = ", concatenate_depths)
    if normalise:
        normalised_y = pd.DataFrame(y) # normalising y
        normalise_standard(normalised_y)
        y = normalised_y.to_numpy().ravel()
        
        # We don't normalise the variance vector
        y_var = np.array(y_var).ravel()

        
    else:
        y = np.array(y).ravel()
        y_var = np.array(y_var).ravel()
    
    
        
    for ker in kernels:
        print(">> Testing GraKeL kernel:", ker)
        mean_sq_error = []
        mean_sq_training_error = []
        for run in range(number_of_runs):
            # Setting up the slices for kfold
            indices = range(len(G))

#                     n_splits = 10
            n_splits = len(G)
            kfold = KFold(n_splits=n_splits, shuffle=True, random_state=run+1)

            concatenated_y_training_pred = ()
            concatenated_y_training_real = ()
            concatenated_y_pred = ()
            concatenated_y_real = ()
            concatenated_indexes = ()
            concatenated_y_pred_variance = ()


            normalise_X = True
            delete_few_occurances = False
            scaler = preprocessing.MinMaxScaler()
#                     scaler = preprocessing.StandardScaler()
            if len(extra_tests) == 0:


                for index_train, index_test in kfold.split(indices):

                    if ker != "RUnits":
                        gk = GRAKEL_KERNELS[ker]() 
                        if not explicit_features:
                            G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=size_of_test, random_state=run)
                            K_train = gk.fit_transform(G_train)
                            K_test = gk.transform(G_test)

                            if len(extra_tests) != 0:
                                K_extra_test = gk.transform(G_extra_tests)
                        else:
                            if printing_results:
                                print(">>> Getting Explicit Vectors from Kernel")
                            param = "temp/dataset_explicit_vectors/" + ker + str(t) + str(len(G)) + "Chiral_" + str(CHIRALITY) + 'final'
                            try:
                                feat = pickle.load(open('{0}.pickle'.format(param), "rb"))
                                if printing_results:
                                    print(">>>> Pickle found")
                            except (OSError, IOError) as e:
                                print(">>>> Pickle not found, recomputing")
                                gk.fit(G)
                                gk.parse_input(G)
                                feat = gk._feat_vectors
                                pickle.dump(feat, open('{0}.pickle'.format(param), "wb"))

#                                     gk.fit(G)
#                                     gk.parse_input(G)
#                                     feat = gk._feat_vectors
                            depth = int(ker[-1]) # this is the last character of the kernel, which indicates depth.
                            if concatenate_depths:
                                feature_vectors = feat[0]
                                for i in range(1, depth +1):
                                    feature_vectors = np.concatenate((feature_vectors, feat[i]), axis=1)
                            else:
                                feature_vectors = feat[depth]
#                                     print(len(gk._inv_labels[depth]))
                            if printing_results:
                                print("Dimension of Feature Vectors: ", len(feature_vectors[0]))
#                                     normalise_X = True

                            if with_shortest_paths:
                                feature_vectors =  np.concatenate((feature_vectors, features_shortest_paths), axis=1)
                            
                            if delete_few_occurances:
#                                 print(feature_vectors)
                                feature_vectors[feature_vectors < 4] = 0
#                                 print(feature_vectors)
                            if normalise_X:
#                                         min_max_scaler = preprocessing.MinMaxScaler()
                                feature_vectors = scaler.fit_transform(feature_vectors)
#                                     # Using PCA to reduce dimensions
#                                     print("Reducing Dimensions")
#                                     feature_vectors = PCA(n_components=20).fit_transform(feature_vectors)


                            if covar_feature_selection:
                                cor = corr2_coeff(feature_vectors.T,feature_vectors.T)
                                p = np.argwhere(np.triu(np.isclose(corr2_coeff(feature_vectors.T,feature_vectors.T),1),1))
                                if printing_results:
                                    print("length before selection:", len(feature_vectors[0]))
                                feature_vectors = np.delete(feature_vectors,p[:,1],axis=1)
                                if printing_results:
                                    print("length after selection:", len(feature_vectors[0]))
                           
                            K_train = feature_vectors[index_train]
                            K_test = feature_vectors[index_test]
                            y_train = y[index_train]
                            y_test = y[index_test]
                            y_var_train = y_var[index_train]

                            if len(extra_tests) != 0:
                                gk.fit(G_extra_tests)
                                gk.parse_input(G_extra_tests)
                                feat_extra = gk._feat_vectors
                                feature_vectors_extra = feat_extra[depth]
                                if normalise_X:
                                    K_extra_test = scaler.transform(feature_vectors_extra)
#                                         K_extra_test = gk.transform(G_extra_tests)



                
                
                    # We define the learning method for each fold. That only makes a difference when using y_var, as we want to make sure we get y_var[index_train] for each fold
                    met = defining_learning_method(t, learning_met, weight_y_var_train, weight_y_var_train_for_tg_and_tg2, y_var[index_train])


                    if RFE_selection:
                        print(">> Performing REF selection for ", RFE_selection_n_features, " features.")
                        selector = RFE(met, n_features_to_select = RFE_selection_n_features, step=0.01)
                        selector = selector.fit(K_train, y_train)
                        feature_vectors = feature_vectors[:, selector.support_]
                        K_train = feature_vectors[index_train]
                        K_test = feature_vectors[index_test]

                    if log_y:
                        met.fit(K_train, np.log(y_train))
                        y_pred, y_pred_variance = met.predict(K_test,  return_std=True)  
                        print("y_pred ", np.exp(y_pred))
                        print("y_test: ", y_test)
                        y_training_pred, y_training_pred_variance = met.predict(K_train,  return_std=True)
                        print("y_training_pred: ", np.exp(y_training_pred))
                        print("y_training_test: ", y_train)
                        print("MSE of train: ", mean_squared_error(y_train, np.exp(y_training_pred)))
                    else:
                        met.fit(K_train, y_train)
                        y_pred, y_pred_variance = met.predict(K_test,  return_std=True)  
                        if truncate_to_zero_if_neg and t in ['E_(MPa)', 'σbreak (MPa)', 'εbreak (pct)']:
                            y_pred[y_pred < 0] = 0
                        if printing_results:
                            print("y_pred ", y_pred)
                            print("y_test: ", y_test)
                        y_training_pred, y_training_pred_variance = met.predict(K_train,  return_std=True)
#                                 print("y_training_pred: ", y_training_pred)
#                                 print("y_training_test: ", y_train)
                        if printing_results:
                            print("MSE of train: ", mean_squared_error(y_train, y_training_pred))
                    if t == 'N_Tgs':
                        print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
#                         if run == 0 :
#                             print("Real:", y_test)
#                             print("Pred:", np.around(y_pred, decimals = 1))
                    if error == 'mean_sq':       
                        if log_y:
                            mean_sq_error.append(mean_squared_error(y_test, np.exp(y_pred)))
                             # also calculating the training error: 
                            mean_sq_training_error.append(mean_squared_error(y_train, np.exp(y_training_pred)))
                        else:
                            mean_sq_error.append(mean_squared_error(y_test, y_pred))
                             # also calculating the training error: 
                            mean_sq_training_error.append(mean_squared_error(y_train, y_training_pred))
                    if error == 'r2':
                        mean_sq_error.append(r2_score(y_test, y_pred))
                        mean_sq_training_error.append(r2_score(y_train, y_training_pred))
                    # for extra tests 
                    if len(extra_tests) != 0:
                        print('Extras')
                        y_extra_pred = met.predict(K_extra_test, return_std=True)
                        for test in range(len(G_extra_tests)):
                            print(block_ratios_to_test[test], y_extra_pred[test])
#                         print(y_extra_pred)



                    #  generating some explanations
#                     graph_to_be_explained = 70
                    graph_to_be_explained = graph_being_explained
#                     graph_to_be_explained = 15
                    if printing_results:
                        print("index_test_is", index_test)
                        print("graph is", numbers_of_polymers[index_test[0]])
                    if explanation and numbers_of_polymers[index_test[0]] == graph_to_be_explained:
                        #  This counts the number of different values each feature can have. If this number is less than the threashold, then 
                        # ... the algorithm will treat them as discrete features, as oppposed to getting quartiles. 
                        threashold = 1

                        # Adding this in case the feature vectors were pickled before. 
                        gk.fit(G)
                        gk.parse_input(G)

                        categorical_features = np.argwhere(np.array([len(set(pd.DataFrame(feature_vectors).iloc[:,x])) for x in range(pd.DataFrame(feature_vectors).shape[1])]) <= threashold).flatten()
                        print("Number of categorical features", len(categorical_features))
#                         print(gk._inv_labels[0]['C'])
#                         print(gk._inv_labels[0]['S'])
#                                 print(gk._inv_labels[depth])
                        list_of_feature_names = list(gk._inv_labels[depth].keys())
                        if RFE_selection:
                            list_of_feature_names = [list_of_feature_names[a] for a in range(len(list_of_feature_names)) if selector.support_[a]]
                        explainer = lime.lime_tabular.LimeTabularExplainer(np.asarray(K_train), 
                                                   feature_names=list_of_feature_names,
                                                   class_names=[ker], 
                                                   categorical_features=categorical_features, 
                                                   verbose=True, 
                                                   mode='regression')

#                         

                        i = 0
#                                 print("K_test is: ", K_test)
                        print("Graph being explained is", names_of_polymers[index_test[i]])
                        mols = []
                        smis = []
#                                 exp = explainer.explain_instance(np.asarray(K_test)[i], met.predict, num_features=5)
                        exp = explainer.explain_instance(K_test[i], met.predict, num_features=50)
                        counter = 0 # to know whether it is the first time we see a motiff on test set
                       
                        print(gk._patterns)
                        for importance in range(50):
                            node, some_test_graph, pattern, found_in_test_set, change_in_proba = generate_explanations_and_subgraphs(exp, gk, importance, i, depth, index_test, index_train)

                            print("node: ", node)
                            print("graph: ", some_test_graph)
                            print("graph's database number: ", numbers_of_polymers[some_test_graph])
                            print("graph's database name: ", names_of_polymers[some_test_graph])
                            # Generating the subgraph assocaited to center node and depth found in a test graph. 
                            nx_subgraph = generate_networkx_subgraph_from_grakel_graph_node_and_depth(G[some_test_graph], depth, node)
                            # Getting its smiles strings
                            smi = pysmiles.write_smiles(nx_subgraph, default_element='*', start=None)
                            print(smi)
                            # Transforming the smiles string in a nice looking chemistry diagram. 
                            mol = Chem.MolFromSmiles(smi)
#                             print(mol)
                            mols.append(mol)
                            smis.append(smi)
                            if found_in_test_set:
                                if counter == 0:
#                                     add_row_to_csv_file_of_smiles([change_in_proba, smi], 's_break_importances.csv', True)
                                    add_row_to_csv_file_of_smiles([change_in_proba, smi], name_file_explanations, True)
                                    counter += 1
                                else:
                                    add_row_to_csv_file_of_smiles([change_in_proba, smi], name_file_explanations, False)
#                                     add_row_to_csv_file_of_smiles([change_in_proba, smi], 's_break_importances.csv', False)
#                             return gk._patterns
#                             return mols, smis

                    concatenated_indexes = list(concatenated_indexes) + list(index_test.ravel())
                    if log_y:
                        concatenated_y_pred = list(concatenated_y_pred) + list(np.exp(y_pred).ravel())
                    else: 
                        concatenated_y_pred = list(concatenated_y_pred) + list(y_pred.ravel())
                    concatenated_y_real = list(concatenated_y_real) + list(y_test.ravel())
            #         print(concatenated_y_pred, concatenated_y_real)

                    concatenated_y_pred_variance = list(concatenated_y_pred_variance) + list(y_pred_variance.ravel())
                    if log_y:
                        concatenated_y_training_pred = list(concatenated_y_training_pred) + list(np.exp(y_training_pred).ravel())
                    else:
                        concatenated_y_training_pred = list(concatenated_y_training_pred) + list(y_training_pred.ravel())

                    concatenated_y_training_real = list(concatenated_y_training_real) + list(y_train.ravel())

                if visualisation and number_of_runs < 2:
#                  
                    if len(extra_tests) == 0: # if there is an indexes vector
                        plot_predictions_against_real_values_with_hovering(concatenated_y_real, concatenated_y_pred, concatenated_y_pred_variance, y_var, ker, t, concatenated_indexes, names_of_polymers, numbers_of_polymers, colours_of_polymers)
                    else:
                        plot_predictions_against_real_values(y_test, y_pred, ker, t)



        
        # Saving the average errros
        print("> Saving Results")
        mean = np.mean(mean_sq_error)
        std = np.std(mean_sq_error)
        print('Predicting', t, 'For', ker, 'MSE is, on average =', round(mean,4), 'and std dev of MSE =', round(std, 4))
        mean_training = np.mean(mean_sq_training_error)
        std_training = np.std(mean_sq_training_error)            

        print('Training', t, 'For', ker, 'MSE is, on average =', round(mean_training,4), 'and std dev of MSE =', round(std_training, 4))
        #  Saving values in a pickled file. 
        parameters = [t, ker, learning_met, str(number_of_runs), str(explicit_features), str(NORMALIZING_GRAPH_KERNELS), str(normalise), str(training_error), str(weight_y_var_train), str(weight_y_var_train_for_tg_and_tg2)]

        # With Training Error 
        if training_error:
            print(">> Saving traning error")
#                 parameters = [t, ker, learning_met, str(number_of_runs), str(explicit_features), str(NORMALIZING_GRAPH_KERNELS), str(normalise), str(training_error), str(weight_y_var_train), str(weight_y_var_train_for_tg_and_tg2)]
#                 t_error = mean_squared_error(y, y_pred_training_error)
            filename = ''.join(parameters)
            path = 'temp/results/' + filename
            pickle.dump([parameters, mean, std, mean_training],open('{0}.pickle'.format(path),'wb'))

        else:
            print(">> Not Saving traning error")
            filename = ''.join(parameters)
            path = 'temp/results/' + filename
            pickle.dump([parameters, mean, std],open('{0}.pickle'.format(path),'wb'))
            
            
param = "temp/dataset_explicit_vectors/" + "spath_feat" + "all_graphs_" + "Chiral_True" 
features_sp = pickle.load(open('{0}.pickle'.format(param), "rb"))            