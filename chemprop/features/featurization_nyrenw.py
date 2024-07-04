from argparse import Namespace
from typing import List, Tuple, Union
from rdkit import Chem
import torch
from torch.utils.data import Dataset

import h5py
import numpy as np
import pandas as pd
import time
import os
# Atom feature sizes
MAX_ATOMIC_NUM = 20
ATOM_FEATURES = {
    'atomic_num': list(range(MAX_ATOMIC_NUM)),
    'degree': [0, 1, 2, 3, 4, 5],
    'formal_charge': [-1, -2, 1, 2, 0],
    'chiral_tag': [0, 1, 2, 3],
    'num_Hs': [0, 1, 2, 3, 4],
    'hybridization': [
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3,
        Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2
    ],
}

# Distance feature sizes
PATH_DISTANCE_BINS = list(range(10))
THREE_D_DISTANCE_MAX = 20
THREE_D_DISTANCE_STEP = 1
THREE_D_DISTANCE_BINS = list(range(0, THREE_D_DISTANCE_MAX + 1, THREE_D_DISTANCE_STEP))

# len(choices) + 1 to include room for uncommon values; + 2 at end for IsAromatic and mass
ATOM_FDIM = sum(len(choices) + 1 for choices in ATOM_FEATURES.values()) + 2
BOND_FDIM = 14

# Memoization
SMILES_TO_GRAPH = {}


def clear_cache():
    """Clears featurization cache."""
    global SMILES_TO_GRAPH
    SMILES_TO_GRAPH = {}


def get_atom_fdim(args: Namespace) -> int:
    """
    Gets the dimensionality of atom features.

    :param: Arguments.
    """
    return ATOM_FDIM


def get_bond_fdim(args: Namespace) -> int:
    """
    Gets the dimensionality of bond features.

    :param: Arguments.
    """
    return BOND_FDIM


def onek_encoding_unk(value: int, choices: List[int]) -> List[int]:
    """
    Creates a one-hot encoding.

    :param value: The value for which the encoding should be one.
    :param choices: A list of possible values.
    :return: A one-hot encoding of the value in a list of length len(choices) + 1.
    If value is not in the list of choices, then the final element in the encoding is 1.
    """
    encoding = [0] * (len(choices) + 1)
    index = choices.index(value) if value in choices else -1
    encoding[index] = 1

    return encoding


def atom_features(atom: Chem.rdchem.Atom, functional_groups: List[int] = None) -> List[Union[bool, int, float]]:
    """
    Builds a feature vector for an atom.

    :param atom: An RDKit atom.
    :param functional_groups: A k-hot vector indicating the functional groups the atom belongs to.
    :return: A list containing the atom features.
    """
    features = onek_encoding_unk(atom.GetAtomicNum() - 1, ATOM_FEATURES['atomic_num']) + \
           onek_encoding_unk(atom.GetTotalDegree(), ATOM_FEATURES['degree']) + \
           onek_encoding_unk(atom.GetFormalCharge(), ATOM_FEATURES['formal_charge']) + \
           onek_encoding_unk(int(atom.GetChiralTag()), ATOM_FEATURES['chiral_tag']) + \
           onek_encoding_unk(int(atom.GetTotalNumHs()), ATOM_FEATURES['num_Hs']) + \
           onek_encoding_unk(int(atom.GetHybridization()), ATOM_FEATURES['hybridization']) + \
           [1 if atom.GetIsAromatic() else 0] + \
           [atom.GetMass() * 0.01]  # scaled to about the same range as other features
    if functional_groups is not None:
        features += functional_groups
    return features


def bond_features(bond: Chem.rdchem.Bond) -> List[Union[bool, int, float]]:
    """
    Builds a feature vector for a bond.

    :param bond: A RDKit bond.
    :return: A list containing the bond features.
    """
    if bond is None:
        fbond = [1] + [0] * (BOND_FDIM - 1)
    else:
        bt = bond.GetBondType()
        fbond = [
            0,  # bond is not None
            bt == Chem.rdchem.BondType.SINGLE,
            bt == Chem.rdchem.BondType.DOUBLE,
            bt == Chem.rdchem.BondType.TRIPLE,
            bt == Chem.rdchem.BondType.AROMATIC,
            (bond.GetIsConjugated() if bt is not None else 0),
            (bond.IsInRing() if bt is not None else 0)
        ]
        fbond += onek_encoding_unk(int(bond.GetStereo()), list(range(6)))
    return fbond


class MolGraph:
    """
    A MolGraph represents the graph structure and featurization of a single molecule.

    A MolGraph computes the following attributes:
    - smiles: Smiles string.
    - n_atoms: The number of atoms in the molecule.
    - n_bonds: The number of bonds in the molecule.
    - f_atoms: A mapping from an atom index to a list atom features.
    - f_bonds: A mapping from a bond index to a list of bond features.
    - a2b: A mapping from an atom index to a list of incoming bond indices.
    - b2a: A mapping from a bond index to the index of the atom the bond originates from.
    - b2revb: A mapping from a bond index to the index of the reverse bond.
    """

    def __init__(self, smiles: str, args: Namespace, arr: np.array = None):
        """
        Computes the graph structure and featurization of a molecule.

        :param smiles: A smiles string.
        :param args: Arguments.
        """
        self.smiles = smiles
        self.target = arr
        self.n_atoms = 0  # number of atoms
        self.n_bonds = 0  # number of bonds
        self.f_atoms = []  # mapping from atom index to atom features
        self.f_bonds = []  # mapping from bond index to concat(in_atom, bond) features
        self.a2b = []  # mapping from atom index to incoming bond indices
        self.b2a = []  # mapping from bond index to the index of the atom the bond is coming from
        self.b2revb = []  # mapping from bond index to the index of the reverse bond
        # Convert smiles to molecule
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol) # MF 020620

        # fake the number of "atoms" if we are collapsing substructures
        self.n_atoms = mol.GetNumAtoms()
        
        # Get atom features
        for i, atom in enumerate(mol.GetAtoms()):
            self.f_atoms.append(atom_features(atom))
        self.f_atoms = [self.f_atoms[i] for i in range(self.n_atoms)]

        for _ in range(self.n_atoms):
            self.a2b.append([])

        # Get bond features
        for a1 in range(self.n_atoms):
            for a2 in range(a1 + 1, self.n_atoms):
                bond = mol.GetBondBetweenAtoms(a1, a2)

                if bond is None:
                    continue

                f_bond = bond_features(bond)

                if args.atom_messages:
                    self.f_bonds.append(f_bond)
                    self.f_bonds.append(f_bond)
                else:
                    self.f_bonds.append(self.f_atoms[a1] + f_bond)
                    self.f_bonds.append(self.f_atoms[a2] + f_bond)

                # Update index mappings
                b1 = self.n_bonds
                b2 = b1 + 1
                self.a2b[a2].append(b1)  # b1 = a1 --> a2
                self.b2a.append(a1)
                self.a2b[a1].append(b2)  # b2 = a2 --> a1
                self.b2a.append(a2)
                self.b2revb.append(b2)
                self.b2revb.append(b1)
                self.n_bonds += 2
                
        self.max_num_bonds = max(0, max(len(in_bonds) for in_bonds in self.a2b))
        self.a2b = [self.a2b[a] + [np.nan] * (self.max_num_bonds - len(self.a2b[a])) for a in range(self.n_atoms)]

        
    
                
    def save_to_h5py(self, f: h5py.File, idx: int):
        """
        Saves the MolGraph to an HDF5 file.

        :param f: An HDF5 file.
        :param idx: The index at which to save the MolGraph.
        """
        if 'num_mols' not in f:
            f.create_dataset('num_mols', data=np.array([0]))
            
        # Save MolGraph variables to h5py file.abs
        # Use three structure with top level nodes of index 0-1000, 1000-2000, 2000-3000...
        idx_tmp = int(idx/1000)
        idx = idx % 1000

        # if group idx_tmp does not exist, create it
        if str(idx_tmp) not in f:
            top_group = f.create_group(str(idx_tmp), track_order=True)
        else:
            top_group = f[str(idx_tmp)]
        # create subgroup for this molecule
        group = top_group.create_group(str(idx), track_order=True)
        group.create_dataset('f_atoms', data=self.f_atoms, compression="lzf")
        group.create_dataset('f_bonds', data=self.f_bonds, compression="lzf")
        group.create_dataset('a2b', data=self.a2b, compression="lzf")
        group.create_dataset('b2a', data=self.b2a, compression="lzf")
        group.create_dataset('b2revb', data=self.b2revb, compression="lzf")
        group.create_dataset('n_atoms_bonds', data=np.array([self.n_atoms, self.n_bonds]), compression="lzf")
        group.create_dataset('smiles', data=np.char.encode(self.smiles, encoding='cp037'))
        group.create_dataset('target', data=self.target, compression="lzf")
        f['num_mols'][0] += 1
        return
    
class BatchMolGraphH5py(Dataset):
    def __init__(self, h5py_file: str, args: Namespace):
        self.h5py_file = h5py.File(h5py_file, 'r', libver='latest', swmr=True)
        self.smiles_batch = []
        self.n_mols_tot = self.h5py_file['num_mols'][0]
        self.n_mols = 0
        self.n_atoms = 1  # number of atoms
        self.n_bonds = 1  # number of bonds
        self.atom_fdim = get_atom_fdim(args)
        self.bond_fdim = get_bond_fdim(args) + (not args.atom_messages) * self.atom_fdim
        
    def __len__(self):
        return self.n_mols_tot
    
    def __getitem__(self, idx):
        """
        idx_tmp = int(idx/1000)
        idx = idx % 1000
        group = self.h5py_file[str(idx_tmp)][str(idx)]
        datasets = group.items()
        datasets = [d for d in datasets]
        f_atoms = np.float32(datasets[0][1][:])
        f_bonds =  np.float32(datasets[1][1][:])
        arr =  np.float32(datasets[2][1][:])
        b2a =  np.float32(datasets[3][1][:])
        b2revb =  np.float32(datasets[4][1][:])
        n_atoms_bonds =  datasets[5][1][:]
        smiles = np.char.decode(datasets[6][1][...], encoding='cp037').astype(str)
        targets =  np.float32(datasets[7][1][:])
        return f_atoms, f_bonds, arr, b2a, b2revb, n_atoms_bonds, smiles, targets
        """
        return idx

    def __getitems__(self, idx_lst):
        n_mols = len(idx_lst)
        self.smiles_batch = []
        self_n_atoms = 1  # number of atoms
        self_n_bonds = 1  # number of bonds
        a_scope = []  # list of tuples indicating (start_atom_index, num_atoms) for each molecule
        b_scope = []  # list of tuples indicating (start_bond_index, num_bonds) for each molecule
        a2b = []  # mapping from atom index to incoming bond indices
        b2a = [0]  # mapping from bond index to the index of the atom the bond is coming from
        b2revb = [0]  # mapping from bond index to the index of the reverse bond
        f_atoms = np.array([[0] * self.atom_fdim])  # atom features
        f_bonds = np.array([[0] * self.bond_fdim])  # combined atom/bond features
        targets = np.zeros((n_mols, 50))
        for i, idx in enumerate(idx_lst):
            idx_tmp = int(idx/1000)
            idx = idx % 1000
            group = self.h5py_file[str(idx_tmp)][str(idx)]
            datasets = group.items()
            datasets = [d for d in datasets]
            f_atoms = np.concatenate((f_atoms, datasets[0][1][:]), axis=0)
            f_bonds = np.concatenate((f_bonds, datasets[1][1][:]), axis=0)
            arr = datasets[2][1][:] + self_n_bonds
            a2b.append(arr)
            b2a = np.concatenate((b2a, datasets[3][1][:] + self_n_atoms), axis=0)
            b2revb = np.concatenate((b2revb, datasets[4][1][:] + self_n_bonds), axis=0)
            n_atoms_bonds = datasets[5][1][:]
            n_atoms = n_atoms_bonds[0]
            n_bonds = n_atoms_bonds[1]
            a_scope.append((self_n_atoms, n_atoms))
            b_scope.append((self_n_bonds, n_bonds))
            smiles = np.char.decode(datasets[6][1][...], encoding='cp037').astype(str)
            self.smiles_batch.append(smiles)
            self_n_atoms += n_atoms
            self_n_bonds += n_bonds
            targets[i][:] = datasets[7][1][:]
        max_num_bonds = max(1, max(in_bonds.shape[1] for in_bonds in a2b)) # max with 1 to fix a crash in rare case of all single-heavy-atom mols
        a2b = np.vstack([np.concatenate((arr, np.zeros((arr.shape[0], max_num_bonds - arr.shape[1]))), axis=1) for arr in a2b])
        # Pad with zeros on the top row of a2b
        a2b = np.nan_to_num(np.concatenate((np.zeros((1, max_num_bonds)), a2b), axis=0))
        
        
        # Convert to torch tensors
        f_atoms = torch.FloatTensor(f_atoms)
        f_bonds = torch.FloatTensor(f_bonds)
        a2b = torch.LongTensor(a2b)
        b2a = torch.LongTensor(b2a)
        b2revb = torch.LongTensor(b2revb)
        targets = torch.FloatTensor(targets)


        return f_atoms, f_bonds, a2b, b2a, b2revb, a_scope, b_scope, targets   
    def __del__(self):
        self.h5py_file.close()
        return
    
    def collate_fn(*input):
        input = input[1]
        return input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7]
    #def collate_fn(self, *input):
    #    print(type(input[0][0][0]))
    #    print(torch.utils.data.get_worker_info())
    
    def get_b2b(self) -> torch.LongTensor:
        """
        Computes (if necessary) and returns a mapping from each bond index to all the incoming bond indices.

        :return: A PyTorch tensor containing the mapping from each bond index to all the incoming bond indices.
        """

        if self.b2b is None:
            b2b = self.a2b[self.b2a]

class BatchMolGraph:
    """
    A BatchMolGraph represents the graph structure and featurization of a batch of molecules.

    A BatchMolGraph contains the attributes of a MolGraph plus:
    - smiles_batch: A list of smiles strings.
    - n_mols: The number of molecules in the batch.
    - atom_fdim: The dimensionality of the atom features.
    - bond_fdim: The dimensionality of the bond features (technically the combined atom/bond features).
    - a_scope: A list of tuples indicating the start and end atom indices for each molecule.
    - b_scope: A list of tuples indicating the start and end bond indices for each molecule.
    - max_num_bonds: The maximum number of bonds neighboring an atom in this batch.
    - b2b: (Optional) A mapping from a bond index to incoming bond indices.
    - a2a: (Optional): A mapping from an atom index to neighboring atom indices.
    """

    def __init__(self, mol_graphs: List[MolGraph], args: Namespace):
        self.smiles_batch = [mol_graph.smiles for mol_graph in mol_graphs]
        self.n_mols = len(self.smiles_batch)

        self.atom_fdim = get_atom_fdim(args)
        self.bond_fdim = get_bond_fdim(args) + (not args.atom_messages) * self.atom_fdim

        # Start n_atoms and n_bonds at 1 b/c zero padding
        self.n_atoms = 1  # number of atoms (start at 1 b/c need index 0 as padding)
        self.n_bonds = 1  # number of bonds (start at 1 b/c need index 0 as padding)
        self.a_scope = []  # list of tuples indicating (start_atom_index, num_atoms) for each molecule
        self.b_scope = []  # list of tuples indicating (start_bond_index, num_bonds) for each molecule

        # All start with zero padding so that indexing with zero padding returns zeros
        f_atoms = [[0] * self.atom_fdim]  # atom features
        f_bonds = [[0] * self.bond_fdim]  # combined atom/bond features
        a2b = [[]]  # mapping from atom index to incoming bond indices
        b2a = [0]  # mapping from bond index to the index of the atom the bond is coming from
        b2revb = [0]  # mapping from bond index to the index of the reverse bond
        for mol_graph in mol_graphs:
            f_atoms.extend(mol_graph.f_atoms)
            f_bonds.extend(mol_graph.f_bonds)

            for a in range(mol_graph.n_atoms):
                a2b.append([b + self.n_bonds for b in mol_graph.a2b[a]])
            for b in range(mol_graph.n_bonds):
                b2a.append(self.n_atoms + mol_graph.b2a[b])
                b2revb.append(self.n_bonds + mol_graph.b2revb[b])
            self.a_scope.append((self.n_atoms, mol_graph.n_atoms))
            self.b_scope.append((self.n_bonds, mol_graph.n_bonds))
            self.n_atoms += mol_graph.n_atoms
            self.n_bonds += mol_graph.n_bonds

        self.max_num_bonds = max(1, max(len(in_bonds) for in_bonds in a2b)) # max with 1 to fix a crash in rare case of all single-heavy-atom mols

        self.f_atoms = torch.FloatTensor(f_atoms)
        self.f_bonds = torch.FloatTensor(f_bonds)
        self.a2b = torch.LongTensor([a2b[a] + [0] * (self.max_num_bonds - len(a2b[a])) for a in range(self.n_atoms)])
        self.b2a = torch.LongTensor(b2a)
        self.b2revb = torch.LongTensor(b2revb)
        self.b2b = None  # try to avoid computing b2b b/c O(n_atoms^3)
        self.a2a = None  # only needed if using atom messages

    def get_components(self) -> Tuple[torch.FloatTensor, torch.FloatTensor,
                                      torch.LongTensor, torch.LongTensor, torch.LongTensor,
                                      List[Tuple[int, int]], List[Tuple[int, int]]]:
        """
        Returns the components of the BatchMolGraph.

        :return: A tuple containing PyTorch tensors with the atom features, bond features, and graph structure
        and two lists indicating the scope of the atoms and bonds (i.e. which molecules they belong to).
        """
        return self.f_atoms, self.f_bonds, self.a2b, self.b2a, self.b2revb, self.a_scope, self.b_scope

    def get_b2b(self) -> torch.LongTensor:
        """
        Computes (if necessary) and returns a mapping from each bond index to all the incoming bond indices.

        :return: A PyTorch tensor containing the mapping from each bond index to all the incoming bond indices.
        """

        if self.b2b is None:
            b2b = self.a2b[self.b2a]  # num_bonds x max_num_bonds
            # b2b includes reverse edge for each bond so need to mask out
            revmask = (b2b != self.b2revb.unsqueeze(1).repeat(1, b2b.size(1))).long()  # num_bonds x max_num_bonds
            self.b2b = b2b * revmask

        return self.b2b

    def get_a2a(self) -> torch.LongTensor:
        """
        Computes (if necessary) and returns a mapping from each atom index to all neighboring atom indices.

        :return: A PyTorch tensor containing the mapping from each bond index to all the incodming bond indices.
        """
        if self.a2a is None:
            # b = a1 --> a2
            # a2b maps a2 to all incoming bonds b
            # b2a maps each bond b to the atom it comes from a1
            # thus b2a[a2b] maps atom a2 to neighboring atoms a1
            self.a2a = self.b2a[self.a2b]  # num_atoms x max_num_bonds

        return self.a2a


def mol2graph(smiles_batch: List[str],
              args: Namespace) -> BatchMolGraph:
    """
    Converts a list of SMILES strings to a BatchMolGraph containing the batch of molecular graphs.

    :param smiles_batch: A list of SMILES strings.
    :param args: Arguments.
    :return: A BatchMolGraph containing the combined molecular graph for the molecules
    """
    mol_graphs = []
    for smiles in smiles_batch:
        if smiles in SMILES_TO_GRAPH:
            mol_graph = SMILES_TO_GRAPH[smiles]
        else:
            mol_graph = MolGraph(smiles, args)
            #if not args.no_cache:
            #    SMILES_TO_GRAPH[smiles] = mol_graph
        mol_graphs.append(mol_graph)
    
    return BatchMolGraph(mol_graphs, args)

def createH5pyFile(df, args: Namespace, h5py_file: str):
    """
    Converts a list of SMILES strings to a BatchMolGraph containing the batch of molecular graphs.

    :param smiles_batch: A list of SMILES strings.
    :param args: Arguments.
    :return: A BatchMolGraph containing the combined molecular graph for the molecules
    """
    # if h5py file exists, skip
    if os.path.exists(h5py_file):
        return
    
    smiles_lst = df['smiles'].tolist()
    target_lst = df['smooth_spectra'].tolist()
    del df
    mol_graphs = []
    f = h5py.File(h5py_file, 'w', libver='latest')
    for ind in range(len(smiles_lst)):
        smiles = smiles_lst.pop(0)
        target = target_lst.pop(0)
        mol_graph = MolGraph(smiles, args, target)
        mol_graph.save_to_h5py(f, ind)
    f.swmr_mode = True
    f.close()
    return
