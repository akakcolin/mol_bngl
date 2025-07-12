#include "define.h"


int **malloc_Mat2d_i(int n1, int n2){
    int i;
    int **array =(int **)malloc(n1*sizeof(int*));
    for(i=0;i<n1;i++){
        array[i]=(int *)malloc(n2*sizeof(int));
    }
    return array;
}

void free_Mat2d_i(int **array, int n1){
  int i;
  for( i=0; i< n1; i++){
    free(array[i]);
  }
  free(array);
}

void copy_ctab_atom(Sdf_Ctab_Atom *a, Sdf_Ctab_Atom *b){
  a->pos[0] = b->pos[0];
  a->pos[1] = b->pos[1];
  a->pos[2] = b->pos[2];
  a->mass_diff = b->mass_diff;
  strcpy(a->symbol, b->symbol);
  a->valence = b->valence;
  a->charge = b->charge;
  a->global_id = b->global_id;
  a->atom_num = b->atom_num;
  a->mol_id = b->mol_id;
}

void copy_ctab_atoms(Sdf_Ctab_Atom *a, Sdf_Ctab_Atom *b, int num_atoms){
  for(int i=0; i<num_atoms; i++){
    if(b[i].global_id > -1){
    copy_ctab_atom(&a[i], &b[i]);
   }
  }
}

void free_molecule(Molecule *mol){
  free(mol->coords);
  free(mol->names);
  free(mol->atom_nums);
  free(mol);
}


double get_dist2(double a[3], double b[3]){
  double v[3];
  v[0] = a[0] - b[0];
  v[1] = a[1] - b[1];
  v[2] = a[2] - b[2];
  return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


int check_site_overlap(Sdf_Ctab_Atom *atoms, int num_atoms, double pos[3], double cutoff_value){
  double r1[3];
  double dist2;
  double cutoff2= cutoff_value*cutoff_value;
  int i;
  Sdf_Ctab_Atom *atom;

   for (int i=0; i< num_atoms; i++){
    atom = atoms+i;
    if(atom->symbol[0]== 'R' && atom->symbol[1] == '#') continue;
    r1[0] = atom->pos[0] - pos[0];
    r1[1] = atom->pos[1] - pos[1];
    r1[2] = atom->pos[2] - pos[2];
    dist2= r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
    if(dist2 < cutoff2){   // too close
      //      printf(" %f %f \n", dist2, cutoff2);
      return 0;
    }
  }
  return 1;
}

int check_mol_overlap(Molecule *mol, Molecule *tmp_mol, double cutoff_scalar) {
  // number of atoms in a molecule
  int mn = mol->num_atoms;
  int cn = tmp_mol->num_atoms;
  int nn = mn + cn;
  double temp;
  int mol_idz = 0;
  double cutoff;
  double rr1[3], rr2[3];
  double dist;
  // Compute inv lattice
  double *radius = malloc(sizeof(double) * nn);

  int i=0;
  int j=0;
  for (i = 0; i < cn; i++) {
    radius[i] = get_pte_vdw(tmp_mol->atom_nums[i]);
  }
  for(i=0; i< mn; i++){
    radius[i+cn] = get_pte_vdw(mol->atom_nums[i]);
  }

   for (i = 0; i < mn; i++) {
    mol_idz = (i + cn);
    rr1[0] = mol->coords[i*3];
    rr1[1] = mol->coords[i*3+1];
    rr1[2] = mol->coords[i*3+2];
    for (j = 0; j < cn; j++) {
      rr2[0] = tmp_mol->coords[j * 3];
      rr2[1] = tmp_mol->coords[j * 3 + 1];
      rr2[2] = tmp_mol->coords[j * 3 + 2];
      cutoff = (radius[j]+radius[mol_idz])*cutoff_scalar;
      cutoff = cutoff * cutoff;
      dist=get_dist2(rr1, rr2);
      if(dist < cutoff){
        //                        printf("dist2 is %f < %f\n", dist, cutoff);
        free(radius);
        return 0;
      }
    }
   }
   free(radius);
   return 1;
}


void next_pos(double *prev, double *pos, double radius, uint64_t *rng){
  // Generate a random point on the sphere
  double dr, u, v, th, phi;
  double dv[3];
  dr = 0.8 + 0.4*rng_uniform(rng);
  dr = dr*radius;
  u = rng_uniform(rng);
  v = rng_uniform(rng);
  th = acos(2.0*v-1.0);
  phi = TWOPI*u;
  dv[0] = dr*sin(th)*cos(phi);
  dv[1] = dr*sin(th)*sin(phi);
  dv[2] = dr*cos(th);
  pos[0] = prev[0] + dv[0];
  pos[1] = prev[1] + dv[1];
  pos[2] = prev[2] + dv[2];
}

void print_mol2xyzfile(Molecule *mol, FILE *outfile) {
  static int counter = 1;

  char *fmt;
  int N = mol->num_atoms;
  fprintf(outfile, " %d\n", N);
  //fprintf(out_file,
  //        "Properties=species:S:1:pos:R:3 pbc=\"F F F\" mol_id=%d "
  //        "struct_number=%d \n",
  //        mol->mol_id, counter);
  fprintf(outfile," Lattice=\"15 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0\" Properties=species:S:1:pos:R:3 pbc=\"T T T\" \n");

  fmt = "%3s % 14.8f % 14.8f % 14.8f\n";
  for (int i = 0; i < N; i++) {
    fprintf(outfile, fmt, atno2sym(mol->atom_nums[i]), mol->coords[i * 3],
            mol->coords[i * 3 + 1], mol->coords[i * 3 + 2]);
  }

  fflush(outfile);
}


Molecule *allocate_molecule(int num_atoms){
  Molecule *mol = malloc(sizeof(Molecule));
  mol->num_atoms = num_atoms;
  mol->names = (int *)malloc(sizeof(int) * num_atoms);
  mol->atom_nums =(int *)malloc(sizeof(int) * num_atoms);
  mol->coords = (double *)malloc(sizeof(double)* num_atoms * 3);
  return mol;
}


Molecule* to_mol(Sdf_Ctab_Atom *atoms, int num_atoms){
  int i;
  Molecule *a=NULL;

  int nl=0;
  Sdf_Ctab_Atom *atom;
  for(i=0; i<num_atoms; i++){
    atom = atoms+i;
    if(atom->symbol[0]== 'R' && atom->symbol[1] == '#') continue;
    nl++;
  }
  a = allocate_molecule(nl);
  a->num_atoms = nl;

  nl = 0;
  for(i=0; i<num_atoms; i++){
    atom = atoms+i;
    a->names[nl] = atom->global_id;
    if(atom->symbol[0]== 'R' && atom->symbol[1] == '#') continue;
    a->names[nl] = atom->global_id ;
    a->atom_nums[nl] =  atom->atom_num;
    memcpy(&a->coords[nl * 3], atom->pos, 3 * sizeof(double));
    nl++;
  }

  return a;
}

// smol is used to append to results

void add_sdf_mol(Sdf_Ctab *smol, Sdf_MetaData *result, Molecule *mol, Settings *set,
                 int sid, int mid){

  Sdf_Ctab *cc=NULL;
  Sdf_Ctab_Counter *counter=NULL;

  int count_inner = 0;
  int verdict = 0;

  //  int fpair, fatom, spair, satom;

  int fatom, satom, l, j, k;
  int ratom1, ratom2;
  int fpair, spair;
  int seed;
  uint64_t rng[2];
  double ranpos[3];
  int num_bonds1, num_bonds2;
  int i;

  int max_total_atoms=0;
  int max_total_bonds=0;
  double pos[3];


  int max_try = set->max_attempts;
  double cutoff_scalar = set->cutoff_scalar;
  double pos_cutoff = set->pos_cutoff;
  double site_cutoff = set->site_cutoff;

  for(i=0; i<result->num_data; i++){
    printf("result current number atoms is\n", result->ctabs[i].counter.num_atoms);
    max_total_atoms += result->ctabs[i].counter.num_atoms;
    max_total_bonds += result->ctabs[i].counter.num_bonds;
  }
  max_total_atoms += smol->counter.num_atoms;
  max_total_bonds += smol->counter.num_bonds;

  srand((unsigned int)time(NULL));
  seed = rand();
  rng_seed(rng, seed);
  rng_next(rng);

  pos[0]=smol->atoms[sid].pos[0];
  pos[1]=smol->atoms[sid].pos[1];
  pos[2]=smol->atoms[sid].pos[2];

  double over1[3];
  int isoverlap=0;
  Molecule *tmp_mol = allocate_molecule(max_total_atoms);

  Sdf_Ctab *tmp_ctab=malloc(sizeof(Sdf_Ctab));
  //Sdf_Ctab *final_ctab = final_meta->ctabs; // just a pointer

  next_pos(pos, over1, pos_cutoff, rng);
  isoverlap = check_site_overlap(smol->atoms, smol->counter.num_atoms, over1, site_cutoff);
  // select one part mol

  sdf_gen_mol(cc2->atoms, cc2->counter.num_atoms, rng);
  over1[0] = over1[0] - cc2->atoms[atom2i - 1].pos[0];
  over1[1] = over1[1] - cc2->atoms[atom2i - 1].pos[1];
  over1[2] = over1[2] - cc2->atoms[atom2i - 1].pos[2];
  for (int ii = 0; ii < cc2->counter.num_atoms; ii++) {
    cc2->atoms[ii].pos[0] += over1[0];
    cc2->atoms[ii].pos[1] += over1[1];
    cc2->atoms[ii].pos[2] += over1[2];
  }
}

typedef struct {
  int num_caps;
  int num_speices_this_type;
  int **cap_index; // [num_speices][cap_index]
  int **unreacted;  // [num_speices][cap_index]
  int cap_type;
} Cap;

void update_mol();
void try_addmol();
int chose_random_mol();
void chose_random_cap();

void accept_rules();

void mc_gen(Sdf_MetaData *meta, BNG_Model *model){
  int num_steps;
  int num_moltypes = model->num_moltypes;
  int i, j, k;
  int *a_or_b=NULL; // {0,1,2,3} vs (different moltypes)
  int *mol_num_species=NULL;
  a_or_b=(int *) malloc(num_moltypes*sizeof(int));
  mol_num_species =(int *)malloc(num_moltypes *sizeof(int));
  for(i=0; i< num_moltypes; i++){
    a_or_b[i]=i;
    mol_num_species[i]=model->species[i].count;
  }

  num_steps = model->species[0].count + model->species[1].count;
  int num_caps=0;
 // initialize cap
  Cap *cap = malloc(sizeof(Cap)* num_moltypes);
  for(i=0; i< model->num_moltypes; i++){
     cap[i].num_caps = model->moltypes[i].num_comp;
     cap[i].cap_index = malloc_Mat2d_i(model->species[i].count, cap[i].num_caps);
     cap[i].unreacted = malloc_Mat2d_i(model->species[i].count, cap[i].num_caps);

     for(k=0; k< model->species[i].count; k++){
       for(j=0; j< cap[i].num_caps; j++){
         cap[i].cap_index[k][j] = model->moltypes[i].components[j].id;
         cap[i].unreacted[k][j] = 1;
       }
     }
     cap[i].cap_type=i; // vs a_or_b
  }

  int chose_mol_type=0;
  int chose_mol=0;
  Sdf_Ctab *lastCtab=NULL;
  // default use first moltype as start
  //  copy_ctab_atoms

  for(i=0; i< num_steps; i++){
    // chose a new mol from species
    int molid=chose_random_mol(model);
    for(int start=0; start < num_moltypes; start++){
      int counter = 0;
      int AB = a_or_b[start];
      int BA=AB+1; // use random generate
      char *a_mol_name=NULL;
      char *b_mol_name=NULL;
      char *a_cap_name=NULL;
      char *b_cap_name=NULL;

      int a_cap_localid, b_cap_localid;
      BNG_BindingRule *rule=NULL;

      if(cap[AB].unreacted[0][0] == -1){
        a_mol_name = model->moltypes[AB].name;
        a_cap_localid=cap[AB].cap_index[0][0]; // component's id
        a_cap_name = model->moltypes[AB].components[0].name;
      }
      if(cap[BA].unreacted[0][0] == -1){
        // random chose a moltype and cap
        b_mol_name = model->moltypes[BA].name;
        b_cap_localid=cap[BA].cap_index[0][0];
        b_cap_name = model->moltypes[BA].components[0].name;   //
      }

      int r=0;
      int cap_allowed=0;
      while(r<model->num_rules){
        rule = &(model->rules[r]);
        if(strcmp(rule->aname, a_cap_name) == 0 && (strcmp(rule->bname, b_cap_name))){
          if(rule->bid[0] == a_cap_localid && rule->bid[1] == b_cap_localid){

            cap_allowed=1;
          }
        }
        r++;
      }
      // into generate new /append to lastCtab
      if(cap_allowed == 1){
        // add mol to lastCtab
        cap[AB].unreacted[0][0] = -1;
        cap[BA].unreacted[0][0] = -1;
      }
    }

    // select one mol's cap
    // check binding rules
    // if allow_binding is True
    // try_to_append to lastMeta
    // update lastMeta
    chose_random_cap();
  }
}

void print_reaction_graph(){

}

void build_network( int *current_network, int *alread_in_network){

}

void cl_primary_external(int *reacted_index, int *cl_index) {

}

void bifunctional_linker( int pos_other){

}
