from rdkit import Chem
import os
import re

import networkx as nx
from itertools import combinations
import matplotlib.pyplot as plt
# parameters part modified base on case

#re_PG_part = re.compile(r'_PG_\(oh~p!*(\%?\d*),oh~s!*(\%?\d*)\)')
#re_TMP_part = re.compile(r'_TMP_\(oh~p!*(\%?\d*),oh~p!*(\%?\d*),oh~p!*(\%?\d*)\)')
#re_ISOPS_part = re.compile(r'_ISOPS_\(cooh!*(\%?\d*),cooh!*(\%?\d*)\)')

re_PG_part = re.compile(r'PG\(oh!*(\%?\d*),oh!*(\%?\d*)\)')
re_TMP_part = re.compile(r'TMP\(oh!*(\%?\d*),oh!*(\%?\d*),oh!*(\%?\d*)\)')
re_ISOPS_part = re.compile(r'ISOPS\(cooh!*(\%?\d*),cooh!*(\%?\d*)\)')



_ISOPS_smile='O=C(O[*])c1cccc(C(=O)O[*])c1'
_PG_smile = 'CC([*])C[*]'
_TMP_smile = 'CCC(C[*])(C[*])C[*]'
_H2O_smile = 'O'

_ISOPS_func_index0=[7,24]
_TMP_func_index0=[7,13,18]
_PG_func_index0=[5,10]

def insert_label_cooh(s, index, label):
    if index <0 or index > len(s):
        return s
    if label=='' or label==None:
        # remove [*]
        return s[:index-2] + s[index+1:]
        #return s
    return s[:index] + ':'+str(label) + s[index:]


def insert_label_oh(s, index, label):
    if index <0 or index > len(s):
        return s
    if label=='' or label==None:
        # remove [*]
        return s[:index-2] + 'O'+ s[index+1:]
        #return s
    return s[:index] + ':'+str(label) + s[index:]

def gen_isops_smile_part(isops, fun_index0, bngl_index):
    # insert from end to begin
    s = isops
    if(fun_index0[0] < fun_index0[1]):
        for index, label in reversed(list(zip(fun_index0, bngl_index))):
            s = insert_label_cooh(s, index, label)
    return s

def gen_tmp_smile_part(isops, fun_index0, bngl_index):
    # insert from end to begin
    s = isops
    if(fun_index0[0] < fun_index0[1]):
        for index, label in reversed(list(zip(fun_index0, bngl_index))):
            s = insert_label_oh(s, index, label)
    return s

def gen_pg_smile_part(isops, fun_index0, bngl_index):
    # insert from end to begin
    s = isops
    if(fun_index0[0] < fun_index0[1]):
        for index, label in reversed(list(zip(fun_index0, bngl_index))):
            s = insert_label_oh(s, index, label)
    return s

def enumerate_one_2dpos(core, *rgroups,randomOrder=False):
    # preserve the positions of the non-dummy core atoms, 
    # we will use these to make sure the cores are drawn
    # the same way in each molecule we generate
    corePos = {}
    conf = core.GetConformer()
    for i in range(conf.GetNumAtoms()):
        corePos[i] = Geometry.Point2D(conf.GetAtomPosition(i))
        
    # Python's itertools handles doing the combinatorics of generating
    # every possible combination of R groups:
    order = itertools.product(*rgroups)
        
    # now we just loop over each combination, copy all the pieces into
    # one molecule, and zip it. That's our product
    for tpl in order:
        tm = Chem.RWMol(core)
        for r in tpl:
            tm.InsertMol(r)
        prod = Chem.molzip(tm)
        if prod is not None:
            # generate 2d coordinates with the core fixed in place
            rdDepictor.Compute2DCoords(prod,canonOrient=False,coordMap=corePos)
            
            # and finally yield the product molecule
            yield prod


def bngl_convert_smiles(bngl_frags):
    ss=None
    frags = bngl_frags.split('.')
    linkers = []
    node_attribs=[]
    color_attr = []
    for idx in range(len(frags)):
        unit_frag = frags[idx]

        if 'ISOPS' in unit_frag:
            res=re.findall(re_ISOPS_part, unit_frag)
            fun_index = _ISOPS_func_index0
            smile_temp = _ISOPS_smile
            new_frag = gen_isops_smile_part(smile_temp, fun_index, res[0])
            color_attr.append('blue')
            node_attribs.append('ISOPS')
        elif 'TMP' in unit_frag:
            res = re.findall(re_TMP_part, unit_frag)
            fun_index = _TMP_func_index0
            smile_temp = _TMP_smile
            new_frag = gen_tmp_smile_part(smile_temp, fun_index, res[0])
            color_attr.append('red')
            node_attribs.append('TMP')
        elif 'PG' in unit_frag:
            res = re.findall(re_PG_part, unit_frag)
            fun_index = _PG_func_index0
            smile_temp = _PG_smile
            new_frag = gen_pg_smile_part(smile_temp, fun_index, res[0])
            color_attr.append('yellow')
            node_attribs.append('PG')
        elif 'H2O' in unit_frag:
            smile_temp = _H2O_smile
            res=['0']
            new_frag =_H2O_smile
            color_attr.append('blue')
            node_attribs.append('H2O')

        else:
            print("not right bnglfrags")
        #new_frag = gen_smile_part(smile_temp, fun_index, res[0])

        linkers.append(res[0])
        print(new_frag)
        if idx ==0:
            ss = new_frag
        else:
            ss = ss + '.'+ new_frag
    return ss, linkers, node_attribs, color_attr


def find_empty_indices(lst, empty_values=(None, '', [], {}, set())):
    return [index for index, element in enumerate(lst) if element in empty_values]

def bngl_convert_graph(linkers, node_attribs):
    total_node = len(linkers)
    G= nx.Graph()
    active_nodes=[]
    for index, link in enumerate(linkers):
        node = index+1
        unpaired = find_empty_indices(link, empty_values=(''))
        if(len(unpaired)>0):
            active_nodes.append([index, unpaired, node_attribs[index]])
        linked_nodes = [n for n in link if n]
        if(linked_nodes):
            #print(node_attribs[index])
            G.add_node(node, name=node_attribs[index])
            for linked_node in linked_nodes:
                for in2, lk2 in enumerate(linkers):
                    if linked_node in lk2 and in2 != index:
                        G.add_edge(node, in2+1)
                        break
    #print(active_nodes)
    return G, active_nodes

def find_angle_topol(G):
    angle_sets = set()
    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        for comb in combinations(neighbors, 2):
            angle_sets.add((node, comb[0], comb[1]))
    return list(angle_sets)


def g2mol(G):
    all=list(G.nodes(data=True))
    atom_symbols = [n[1]['name'] for n in all ]
    atom_index = [n[0] for n in all]
    bonds = [edge for edge in G.edges]
    angle_topol = find_angle_topol(G)

    atom_numbers=[]
    for atom in atom_symbols:
        if atom == 'TMP':
            atom_numbers.append(1)
        else:
            atom_numbers.append(2)

    print(bonds)

    print(angle_topol)

    return {
        'atom_symbols': atom_symbols,
        'atom_types': atom_numbers,
        "bonds": bonds,
        "angles": angle_topol
    }

def write_lmp_topol(doc_title, topol, coords, box):
    with open(doc_title+'.data', 'w') as f:
        header = 'fene2lmp:'
        f.write(f'{header[-220:len(header)]}\n\n')
        f.write('{} atoms\n'.format(len(atoms)))
        f.write('{} bonds\n'.format(len(bonds)))
        f.write('{} angles\n'.format(len(angles)))
        f.write('0 dihedrals\n')
        f.write('0 impropers\n')
        f.write('\n')
        # write structure quantity types
        f.write('{} atom types\n')
        f.write('1 bond types\n')
        f.write('1 angle types\n')
        f.write('0 dihedral types')
        f.write('0 impropers types')

        # write box size
        f.write('{:>12.9f} {:^9.9f} {} {}\n'.format(box.xlo, box.xhi, 'xlo', 'xhi'))
        f.write('{:>12.9f} {:^9.9f} {} {}\n'.format(box.ylo, box.yhi, 'ylo', 'yhi'))
        f.write('{:>12.9f} {:^9.9f} {} {}\n'.format(box.zlo, box.zhi, 'zlo', 'zhi'))
        if box.xy != 0 or box.xz != 0 or box.yz != 0:
            f.write('{:>12.9f} {:^9.9f} {:^9.9f} {} {} {}\n'.format(box.xy, box.xz, box.yz, 'xy', 'xz', 'yz'))


    f.write('\n')
    f.write('\nMasses\n\n')
    for i in masses:
        mass = masses[i]
        ID = '{t:<{s}}'.format(t=str(i), s=3)
        comment = '{:^2} {:5}'.format('#', mass.type)
        f.write('{} {:^12.8f} {:^2}\n'.format(ID, mass.coeffs, comment))


    f.write("\nAtoms\n\n")
    for i in atoms:
        atom = atoms[i]
        atomtype ='{t:<{s}}'.format(t=str(atom_types[i]), s=3)
        comment = '{:^2} {}'.format('#', atom_symbols[i])
        f.write('{:^6} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atomtype, coords[i][0], coords[i][1], coords[i][2], 0, 0, 0, comment))

    f.write("\nBonds\n\n")
    for i in bonds:
        bond = bonds[i]
        bond_type = 1
        bondtype = '{t:<{s}}'.format(t=str(bond_type), s=3)
        id1, id2 = bond
        f.write('{:^5} {} {:^5} {:^5}\n'.format(i, bondtype, id1, id2))

        # Write angles

    f.write("\nAngles\n\n")

    for i in angles:
        angle = angles[i]
        angle_type = 1
        angletype = '{t:<{s}}'.format(t=str(angle_type), s=3)
        id1, id2, id3 = angle
        f.write('{:^5} {} {:^5} {:^5} {:^5}\n'.format(i, angletype, id1, id2, id3))


def sum_active_site(active_nodes):
    sum_oh = 0
    sum_cooh=0
    for node in active_nodes:
        if(node[2] == 'TMP'):
            sum_oh += len(node[1])
        if(node[2] == 'ISOPS'):
            sum_cooh += len(node[1])
    return sum_oh, sum_cooh



def get_polysmile(bngl_ss):
    mol = Chem.MolFromSmiles(bngl_ss)
    mol_com = Chem.molzip(mol)
    newsmi = Chem.MolToSmiles(mol_com)
    #newsmi = Chem.CanonSmiles(newsmi)
    return newsmi


all_species=None
#with open('./example1_reference.species') as f:
with open('./test_10.species') as f:
    l = f.readlines()
    for i in range(len(l)):
        l[i] = l[i].split()
        if l[i][0][0] == '#':
            l[i] = []
    for i in range(len(l)-1,-1,-1):
        if l[i] == []:
            l.pop(i)
    all_species = l

    for i in range(len(l)):
        specie=l[i][0]
        ss,  edges, node_attribs, color_attr = bngl_convert_smiles(specie)
        #print(edges)
        #print(specie)
        G, unpaired_node=bngl_convert_graph(edges, node_attribs)
        topol=g2mol(G)
        oh,cooh = sum_active_site(unpaired_node)
        print('oh:', oh)
        print('cooh:', cooh)
        try:
            cycle = nx.find_cycle(G)
        except:
            cycle = []
        print(cycle)
        #nx.draw(G, node_color=color_attr, with_labels=True)
        #plt.show()
        #name = nx.get_node_attributes(G, 'color')
        #print(name)
        #if(cycle):
        smi = get_polysmile(ss)
        print(smi)

#for specie in all_species:
    #print(specie[0])
#    res=re.findall(re_ISOPS_part, specie[0])
    #print(res)

#def change_string(s, pattern, labels):
#    counter =0
#    def replace(match):
#        nonlocal counter
#        label = labels[counter] if counter < len(labels) else ''
#        counter +=1
#        return pattern + label
#    new_s = re.sub(pattern, replace, s)
#    return new_s


#specie='_ISOPS_(cooh!1,cooh)._TMP_(oh~p!2,oh~p!1,oh~p!3)._ISOPS_(cooh,cooh!2)._ISOPS_(cooh,cooh!3)'
#specie='_ISOPS_(cooh!1,cooh)._PG_(oh~p!1,oh~s!2)._ISOPS_(cooh!2,cooh)'
#res = re.findall(re_TMP_part, specie)
#print(res)

#frags = specie.split('.')
#print("BNGL-species: ", frags)
