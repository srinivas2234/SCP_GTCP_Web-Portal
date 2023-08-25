import sys
import pickle
import shutil
import os
file_name=sys.argv[1]+".txt"
from rdkit import Chem
from rdkit.Chem import Draw
periodic_table=Chem.GetPeriodicTable()
s=sys.argv[1]
print(s)
s=s.split('/')
s=s[2]
destiny_folder="./gtcp_smiles_images_"+s
if(os.path.exists(destiny_folder)):
    shutil.rmtree(destiny_folder)
os.mkdir(destiny_folder)
def MolFromGraphs(node_list,adjacency_matrix,image_name):
        # create empty editable mol object
        mol = Chem.RWMol(Chem.MolFromSmiles(''))
        # add atoms to mol and keep track of index
        node_to_idx = {}
        #print(node_list)
        for i in range(len(node_list)):
            a = Chem.Atom(periodic_table.GetAtomicNumber(str(mapping[str(node_list[i])])))
            #print(type(a))
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

        # Convert RWMol to Mol object
        destination=destiny_folder+'/'+image_name+".png"
        mol = mol.GetMol()
        Draw.MolToFile(mol,destiny_folder+'/'+str(image_name)+".svg")
        
        return mol
fs=open("./graphdata/mapping.txt",'r')
mapping={}
for i in fs:
    r=i.strip('\n')
    r=r.split(" ")
    mapping[r[1]]=r[0]
f=open(file_name,'r')
list_all_graphs=[]
t={}
nodes=[]
edges=[]
for i in f:
    i=i.strip('\n')
    i=i.strip('\r')
    i=i.split(" ")
    if(i[0]=='t'):
        #print("t")
        if(len(t)!=0):
            t["string"]=Chem.MolToSmiles(MolFromGraphs(nodes,edges,t["graph_id"]))
            list_all_graphs.append(t)

        t={}
        t={"graph_id":i[2],"graph":{"nodes":[],"links":[]}}
        nodes=[]
        node_ids=[]
        edges=[]
    elif(i[0]=='v'):
        #print("vertex")
        t["graph"]["nodes"].append({"name":str(i[1]),"value":mapping[str(i[2])]})
        nodes.append(i[2])
        node_ids.append(i[1])
    else:
        #print("edge")
        if(t["graph"]["links"]==[]):
            for j in range(len(nodes)):
                te=[]
                for k in range(len(nodes)):
                    te.append(0)
                edges.append(te)
        t["graph"]["links"].append({"source":str(i[1]),"target":str(i[2])})
        #print("hello")
        edges[node_ids.index(i[1])][node_ids.index(i[2])]=int(i[3])
if(len(t)!=0):
    t["string"]=Chem.MolToSmiles(MolFromGraphs(nodes,edges,t["graph_id"]))
    list_all_graphs.append(t)
print("data_images",len(list_all_graphs))
with open("./all_graphs.txt",'wb') as all:
    pickle.dump(list_all_graphs,all)
