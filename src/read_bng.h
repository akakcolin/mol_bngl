#ifndef READ_BNG_H
#define READ_BNG_H


#include "define.h"

int read_bngl(char *filename, BNG_Model *model);

int parse_molecule_types(char *sources, BNG_MoleculeType *mol);

void print_model(BNG_Model *model);

void free_BNG_Model(BNG_Model *model);

void alloc_BNG_Model(BNG_Model *model, int num_mol, int num_species, int num_rules, int num_bindings);

int parse_species(char *sources, BNG_Species *specie);

#endif
