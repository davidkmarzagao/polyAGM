from pysmiles import read_smiles
import networkx as nx
import numpy as np
import random
import networkx as nx
from rdkit import Chem
import matplotlib.pyplot as plt



def separate_monomers(smile):
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
            block_a = block_a + "}"
            block_b = "{" + block_b
            # print(random_suffling)
            for tpe in random_suffling:
                if tpe == 'a':
                    full_polymer += block_a
                else: full_polymer += block_b
        else:
            for counter in range(round(float(block_ratio[block]))):
                full_polymer += s[block]
    return full_polymer

def translate(s):
    # s is the input Smiles string
    # returns a networkx network, including R/S stereo chemical information
    G = nx.Graph()
    mol = Chem.MolFromSmiles(s)
    id = 0
    for at in mol.GetAtoms():
        G.add_node(id, element = at.GetSymbol())
        id += 1
    
    for bd in mol.GetBonds():
        u, v, w = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx(), bd.GetBondTypeAsDouble()
        G.add_edge(u, v, order = w)
        
    for (id, type) in Chem.FindMolChiralCenters(mol):
        G.nodes[id]['element'] += type
    return G


def from_smiles_to_networkx(smiles, block_ratios, ignore_entries, chiral = False):
    full_pols = []
    for counter in range(len(smiles)):
        if len(ignore_entries) != 0:
            if ignore_entries[counter]:
                continue
        print(counter)
        s = separate_monomers(smiles[counter])
        block_ratio = separate_block_ratio(block_ratios[counter])
        full_polymer = assemble_full_polymer(s, block_ratio)
        if chiral:
            full_pols.append(translate(full_polymer))  # Transforming smiles into networkX
        else:
            # We need to set reinterpret_aromatic to false in order to have double and single bonds alternating in the aromatic ring. 
            # Otherwise all in the ring are set to 1.5 and the WLB algorithm will consider them all single bonds. 
            reint = False
            print("we are considering reinterpret_aromatic", reint)
            full_pols.append(read_smiles(full_polymer, reinterpret_aromatic = reint))  # Transforming smiles into networkX
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
