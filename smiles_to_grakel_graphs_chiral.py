from pysmiles import read_smiles
import networkx as nx
# from .translate import translate
import numpy as np
from rdkit import Chem
import random
import pickle


from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def separate_monomers(smile):
    # print("first input = {}".format(smile))
    split = []
    block = ''
    for char in smile:
        block += char
        if char == '}':
            split.append(block)
            block = ''
    return split

def separate_block_ratio(blocks):
    split = []
    block = ''
    for char in blocks:
        if char == ':':
            split.append(block)
            block = ''
        else:
            block += char
    split.append(block)  # For the last number
    return split


def random_sampling(n_a,n_b):
    # Receive number of monomers of type a, and the number of type b. Returns a list with a shuffle of n_a + n_b items.
    lst = ['a'] * n_a + ['b'] * n_b
    random.shuffle(lst)
    return lst

def assemble_full_polymer(s, block_ratio):
    full_polymer = ''
    for block in range(len(s)):
        if '/' in block_ratio[block]:  # In case we have random sampling of monomers in the middle block
            na, nb = block_ratio[block].split("/")  # the amount of each monomer
            random_suffling = random_sampling(int(na),int(nb))
            block_a, block_b = s[block].split(", ")
            block_a += '}'
            block_b += '{'
            # print("block_a = {}".format(block_a))
            # print("block_b = {}".format(block_b))
                  
            for tpe in random_suffling:
                if tpe == 'a':
                    full_polymer += block_a
                else: full_polymer += block_b
        else:
            for counter in range(round(float(block_ratio[block]))):
                full_polymer += s[block]
    return full_polymer

def from_smiles_to_mol(smiles, block_ratios, ignore_entries):
    print("IN")
    full_mols = []
    for counter in range(len(smiles)):
        if len(ignore_entries) != 0:
            if ignore_entries[counter]:
                continue
        print(counter)
        s = separate_monomers(smiles[counter])
        block_ratio = separate_block_ratio(block_ratios[counter])
        full_polymer = assemble_full_polymer(s, block_ratio)
        full_polymer = full_polymer.replace("{", "")
        full_polymer = full_polymer.replace("}", "")
        full_mols.append(Chem.MolFromSmiles(full_polymer))
    return full_mols

def from_smiles_to_networkx(smiles, block_ratios, ignore_entries, contains_stereo, use_chiral, use_pickle=True, explicit_H = False):
    print("USE CHIRAL = {}".format(use_chiral))
    full_pols = []
    for counter in range(len(smiles)):
        if len(ignore_entries) != 0:
            if ignore_entries[counter]:
                continue
        print(counter)
        s = separate_monomers(smiles[counter])
        block_ratio = separate_block_ratio(block_ratios[counter])
        full_polymer = assemble_full_polymer(s, block_ratio)
        # print('full polymer = {}'.format(full_polymer))
        
        if not use_chiral:
            print("Aromatic false") 
            full_pols.append(read_smiles(full_polymer, reinterpret_aromatic = False))
        else:
            if contains_stereo[counter]:
                print("has stereo")
                if use_pickle:
                    try:
                        cur_pol = pickle.load(open("./temp/pol_networkx/pol_row{}.pickle".format(counter), "rb"))
                        print("loaded")
                    except (OSError, IOError) as e:
                        print("pickle not found, translating")
                        cur_pol = cur_pol = translate(full_polymer)
                        print("dumping")
                        pickle.dump(cur_pol, open("./temp/pol_networkx/pol_row{}.pickle".format(counter), "wb"));
                else: cur_pol = translate(full_polymer)
                full_pols.append(cur_pol)
            else: 
                print("no stereo")
                print("Aromatic false") 
                print("Explicit H", explicit_H) 
                full_pols.append(read_smiles(full_polymer, reinterpret_aromatic = False, explicit_hydrogen = explicit_H))
        # Transforming smiles into networkX
    return full_pols


def from_smiles_to_morgan_dict(smiles, block_ratios, ignore_entries, contains_stereo, use_chiral, use_pickle=True, explicit_H = False, rad = 3):
    print("USE CHIRAL = {}".format(use_chiral))
    full_pols = []
    for counter in range(len(smiles)):
        if len(ignore_entries) != 0:
            if ignore_entries[counter]:
                continue
        print(counter)
        s = separate_monomers(smiles[counter])
        block_ratio = separate_block_ratio(block_ratios[counter])
        full_polymer = assemble_full_polymer(s, block_ratio)
#         print('full polymer = {}'.format(full_polymer))
        full_polymer = full_polymer.replace("{", "")
        full_polymer = full_polymer.replace("}", "")
        if use_pickle:
            try:
                fingerprint = pickle.load(open("./temp/pol_morgan/pol_rad{}_row{}.pickle".format(rad,counter), "rb"))
                print("loaded")
            except (OSError, IOError) as e:
                print("pickle not found, translating")
                mol = Chem.MolFromSmiles(full_polymer)
                generator = rdFingerprintGenerator.GetMorganGenerator(radius = rad, includeChirality=True, fpSize=100000)
                fingerprint = generator.GetSparseCountFingerprint(mol)
#                fingerprint = generator.GetSparCountFingerprint(mol)
                print("dumping")
                pickle.dump(fingerprint, open("./temp/pol_morgan/pol_rad{}_row{}.pickle".format(rad,counter), "wb"));
        else: 
            generator = rdFingerprintGenerator.GetMorganGenerator(radius = rad, includeChirality=True, fpSize=100000)
            fingerprint = generator.GetSparseCountFingerprint(mol)
        full_pols.append(fingerprint)
#
#         mol = Chem.MolFromSmiles(full_polymer)
#         print('full polymer = {}'.format(full_polymer))
#         mol = Chem.MolFromSmiles(full_polymer)
#         if not use_chiral:
#             print("Aromatic false")
#             generator = rdFingerprintGenerator.GetMorganGenerator(radius = rad, includeChirality=True, fpSize=100000)
#             fingerprint = generator.GetSparseCountFingerprint(mol)
# #             fingerprint = generator.GetSparCountFingerprint(mol)
# #             non_zero = fingerprint.GetNonzeroElements()
#             fingerprint
# # non_zero
#             full_pols.append(fingerprint)
#         else:
#             generator = rdFingerprintGenerator.GetMorganGenerator(radius = rad, includeChirality=False, fpSize=100000)
#             print("sparse fingerprint")
#             fingerprint = generator.GetSparseCountFingerprint(mol)
            
#             fingerprint = generator.GetCountFingerprint(mol)
#             non_zero = fingerprint.GetNonzeroElements()        
#          
#             if contains_stereo[counter]:
#                 print("has stereo")
#                 if use_pickle:
#                     try:
#                         cur_pol = pickle.load(open("./temp/pol_networkx/pol_row{}.pickle".format(counter), "rb"))
#                         print("loaded")
#                     except (OSError, IOError) as e:
#                         print("pickle not found, translating")
#                         cur_pol = cur_pol = translate(full_polymer)
#                         print("dumping")
#                         pickle.dump(cur_pol, open("./temp/pol_networkx/pol_row{}.pickle".format(counter), "wb"));
#                 else: cur_pol = translate(full_polymer)
#                 full_pols.append(cur_pol)
#             else: 
#                 print("no stereo")
#                 print("Aromatic false") 
#                 print("Explicit H", explicit_H) 
#                 full_pols.append(read_smiles(full_polymer, reinterpret_aromatic = False, explicit_hydrogen = explicit_H))
        # Transforming smiles into networkX
    return full_pols



def generate_dictionary_of_repeat_units(smiles):
    repeat_uni_dict = dict()
    entry = 0
    for s in smiles:
        for unit in separate_monomers(s):
            if "," in unit:
                for block in unit.split(", "):  # Getting both blocks and then checking if they are alreaddy in the dictionary
                    if "{" not in block:
                        block = "{" + block
                    if "}" not in block:
                        block = block + "}"
                    if block not in repeat_uni_dict:
                        repeat_uni_dict[block] = int(entry)
                        entry += 1
            else:
                if unit not in repeat_uni_dict:
                    repeat_uni_dict[unit] = int(entry)
                    entry += 1
    return repeat_uni_dict

def count_repeat_units(s, block_ratio, repeat_uni_dict):
    counter_repeat_units = np.zeros(len(repeat_uni_dict))
    for block in range(len(s)):
        if '/' in block_ratio[block]: # In case we have random sampling of monomers in the middle block
            na, nb = block_ratio[block].split("/") # the amount of each monomer
            block_a, block_b = s[block].split(", ")
            block_a = block_a + "}"
            block_b = "{" + block_b
            counter_repeat_units[repeat_uni_dict[block_a]] = na
            counter_repeat_units[repeat_uni_dict[block_b]] = nb
        else:
#             print(s[block])
#             print(repeat_uni_dict[s[block]])
            counter_repeat_units[repeat_uni_dict[s[block]]] += round(float(block_ratio[block]))
    return counter_repeat_units

def from_smiles_to_repeat_units(smiles, block_ratios, ignore_entries, repeat_unit_dict):
    full_pols = []
    for counter in range(len(smiles)):
        if ignore_entries[counter]:
            continue
#         print(counter)
        s = separate_monomers(smiles[counter])
        block_ratio = separate_block_ratio(block_ratios[counter])
        full_pols.append(count_repeat_units(s, block_ratio, repeat_unit_dict)) #  Transforming smiles into networkX
    return full_pols


def translate(s):
    # s is the input Smiles string
    # returns a networkx network, including R/S stereo chemical information
    tms = []
    # tms.append(time.process_time())
    ss = ''
    for ch in s:
        if ch != '{' and ch != '}':
            ss += ch;
    # print("ss = {}".format(ss))
    s = ss
    G = nx.Graph()
    # tms.append(time.process_time())
    mol = Chem.MolFromSmiles(s)
    # tms.append(time.process_time())
    id = 0
    for at in mol.GetAtoms():
        G.add_node(id, element=at.GetSymbol())
        id += 1
    # tms.append(time.process_time())
    for bd in mol.GetBonds():
        u, v, w = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx(), bd.GetBondTypeAsDouble()
        G.add_edge(u, v, order=w)
    # tms.append(time.process_time())
    for (id, type) in Chem.FindMolChiralCenters(mol):
        G.nodes[id]['element'] += type    
    # tms.append(time.process_time())
    # print(tms)
    return G