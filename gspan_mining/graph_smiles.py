import os
from rdkit import Chem

#print(Chem.RWMol())
periodic_table=Chem.GetPeriodicTable()
class Smiles():
    def __init__(self):
        self.mapping={}
        f=open("./graphdata/mapping.txt",'r')
        for i in f:
            r=i.strip('\n')
            r=i.split(" ")
            self.mapping[r[1]]=r[0]
    def MolFromGraphs_1(node_list,adjacency_matrix):
        # create empty editable mol object
        mol = Chem.RWMol()

        # add atoms to mol and keep track of index
        node_to_idx = {}
        for i in range(len(node_list)):
            a = Chem.Atom(periodic_table.GetAtomicNumber(str(node_list[i])))
            molIdx = mol.AddAtom(a)
            node_to_idx[i] = molIdx

        # add bonds between adjacent atoms
        for ix, row in enumerate(adjacency_matrix):
            for iy, bond in enumerate(row):

                # only traverse half the matrix
                if iy <= ix:
                    continue

                # add relevant bond type (there are many more of these)
                if bond == 0:
                    continue
                elif bond == 1:
                    bond_type = Chem.rdchem.BondType.SINGLE
                    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
                elif bond == 2:
                    bond_type = Chem.rdchem.BondType.DOUBLE
                    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
                elif bond == 3:
                    bond_type = Chem.rdchem.BondType.TRIPLE
                    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)

        # Convert RWMol to Mol object
        mol = mol.GetMol()            

        return mol
    def MolFromGraphs(self,node_list, adjacency_matrix):
        #print(adjacency_matrix)
        nodes=[]
        for i in node_list:
            nodes.append(self.mapping[i])
        c=Chem.MolToSmiles(self.MolFromGraphs_1(nodes, adjacency_matrix))
        return c