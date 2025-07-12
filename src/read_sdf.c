#include "define.h"
#include "read_bng.h"

// temp method
int get_file_length(const char *filename) {
  FILE *fp;
  char buf[1000];
  int len = 0;
  fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "unable to open %s\n", filename);
    return 0;
  }
  while (fgets(buf, 1000, fp)) {
    len++;
  }
  fclose(fp);
  return len;
}


void concatenateStrings(const char* str1, const char* str2, char* result) {

    int len1 = strlen(str1);
    int len2 = strlen(str2);

    strcpy(result, str1);
    result[len1] = '.';
    strcpy(result + len1 + 1, str2);
}

void read_sdf_v2000_single(Sdf_MetaData *meta, BNG_Model *model) {
  int i;
  FILE *fp =NULL;
  char buffer[LINE_SIZE + 1];
  char *ptr;
  char end_str[4];
  int tmp;
  int records_read=0;

  Sdf_Ctab *cc=NULL;
  Sdf_Ctab_Counter *counter=NULL;

  int num_moltypes = model->num_moltypes;

  meta->num_data = num_moltypes;

  meta->ctabs = malloc(num_moltypes*sizeof(Sdf_Ctab));
  meta->counts = malloc(num_moltypes*sizeof(Sdf_Ctab_Counter));

  for(i=0; i < num_moltypes; i++){
    char filename[100];
    concatenateStrings(model->moltypes[i].name, "sdf", filename);
    fp = fopen(filename, "r");
    if(!fp){
      printf("***ERROR: error in %s file\n", filename);
      exit(EXIT_FAILURE);
    }

    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");

    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");

    cc = meta->ctabs+i;
    cc->num_rgroup=0;

    counter=meta->counts+i;

    records_read =sscanf(buffer, "%d %d %d %d %d %d %d %d %d %d %d %s\n",
                         &(counter->num_atoms), &(counter->num_bonds),
               &(counter->num_atom_lists), &tmp, /* for obsolete record */
               &(counter->flag_chiral), &(counter->num_stext_entries),
               &(counter->num_react_components), &(counter->num_reactants),
               &(counter->num_products), &(counter->num_intermediates),
               &(counter->num_additional_lines), counter->version);

    if (records_read != 10)
      error_exit("Error parsing sdf headline\n");

    counter->last_globalid = -1; // init
    cc->counter.num_atoms = counter->num_atoms;
    cc->counter.num_bonds = counter->num_bonds;

    int num_atoms = counter->num_atoms;
    int num_bonds = counter->num_bonds;
    printf(" index molec %d num_atoms is %d num_bonds is %d\n", i, num_atoms, num_bonds);

    cc->atoms = (Sdf_Ctab_Atom *)malloc(num_atoms*sizeof(Sdf_Ctab_Atom));
    cc->bonds = (Sdf_Ctab_Bond *)malloc(num_bonds*sizeof(Sdf_Ctab_Bond));

    int index_rg=0;
    Sdf_Ctab_Atom *atom=NULL;

    for (int n = 0; n < num_atoms; n++) {
      if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");
      atom = cc->atoms + n;
      ptr = buffer;
      while (isspace(*ptr)) ptr++;
      if (isdigit(*ptr)) {
        records_read =
            sscanf(ptr, "%lf %lf %lf %s %d %d %d %d %d %d %d %d %d %d %d %d",
                   &(atom->pos[0]), &(atom->pos[1]), &(atom->pos[2]), atom->symbol,
                   &(atom->mass_diff), &(atom->charge), &(atom->stereo_parity),
                   &(atom->hydrogen_count), &(atom->flag_stereo_care),
                   &(atom->valence), &(atom->h0_designator),
                   &(atom->type_react_component), &(atom->num_react_component),
                   &(atom->mapping_number), &(atom->flag_invertion),
                   &(atom->flag_exact_change));
        atom->global_id = -1;
        atom->atom_num = atsym2no(atom->symbol); //
        atom->mol_id = i; // biaoji

        if (records_read != 16) error_exit("Error parsing sdf ctab atom data");

        if (atom->symbol[0] == 'R' && atom->symbol[1] == '#') {
          cc->num_rgroup++;
        }
      }
    }
    cc->name_rgroup =(int *)malloc(cc->num_rgroup*sizeof(int));
    cc->atomid_rgroup =(int *)malloc(cc->num_rgroup*sizeof(int));


    Sdf_Ctab_Bond *bond=NULL;

    for(int n=0; n< num_bonds; n++){
      bond = cc->bonds+n;
      if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");
      ptr = buffer;

      while(isspace(*ptr)) ptr++;
      if(isdigit(*ptr)){
        records_read = sscanf(ptr, "%d %d %d %d %d %d",
                              &(bond->fatom),
                              &(bond->satom),
                              &(bond->type),
                              &(bond->stereo),
                              &(bond->topology),
                              &(bond->reacting_center_status)
                              );
        if(records_read != 6) error_exit("Error parsing bond data");
      }
    }

    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");

    // num rgroup in molecule
    int tmp_rgp[9] = {0,0,0,0,0,0,0,0,0};
    if(cc->num_rgroup == 1){
      records_read = sscanf(buffer, "M  RGP %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]));
      if(records_read != 3) error_exit("Error parsing RGP number 1");

      if(tmp_rgp[0] != cc->num_rgroup) error_exit("Error RGroup number in num_rgp=1 \n");
      cc->atomid_rgroup[0] = tmp_rgp[1];  // R's local index
      cc->name_rgroup[0] = tmp_rgp[2]; // R's names
      printf("RGP is %d %d %d\n", cc->num_rgroup,cc->atomid_rgroup[0], cc->name_rgroup[0]);
    }else if(cc->num_rgroup == 2){
      records_read = sscanf(buffer, "M  RGP %d %d %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]),
                            &(tmp_rgp[3]),
                            &(tmp_rgp[4]));
      if(records_read != 5) error_exit("Error parsing RGP number 2");

      if(tmp_rgp[0] != cc->num_rgroup) error_exit("Error RGroup number in num_rgp=2 \n");
      cc->atomid_rgroup[0] = tmp_rgp[1];  // R's local index
      cc->name_rgroup[0] = tmp_rgp[2]; // R's names
      cc->atomid_rgroup[1] = tmp_rgp[3];  // R's local index
      cc->name_rgroup[1] = tmp_rgp[4]; // R's names
      printf("RGP is %d %d %d %d %d \n", cc->num_rgroup, cc->atomid_rgroup[0],
             cc->name_rgroup[0], cc->atomid_rgroup[1], cc->name_rgroup[1]);
    }
    else if(cc->num_rgroup == 3){
       records_read = sscanf(buffer, "M  RGP %d %d %d %d %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]),
                            &(tmp_rgp[3]),
                             &(tmp_rgp[4]),
                             &(tmp_rgp[5]),
                             &(tmp_rgp[6]));
      if(records_read != 7) error_exit("Error parsing RGP number 3");
      if(tmp_rgp[0] != cc->num_rgroup) error_exit("Error RGroup number in num_rgp=3  \n");
      cc->name_rgroup[0] = tmp_rgp[2]; // R's names
      cc->name_rgroup[1] = tmp_rgp[4];
      cc->name_rgroup[2] = tmp_rgp[6];

      cc->atomid_rgroup[0] = tmp_rgp[1]; // R's local index
      cc->atomid_rgroup[1] = tmp_rgp[3];
      cc->atomid_rgroup[2] = tmp_rgp[5];
      printf("RGP is %d %d %d %d %d %d %d \n", cc->num_rgroup,cc->atomid_rgroup[0], cc->name_rgroup[0],
             cc->atomid_rgroup[1], cc->name_rgroup[1], cc->atomid_rgroup[2], cc->name_rgroup[2]);
    }

    else if(cc->num_rgroup == 4){
       records_read = sscanf(buffer, "M  RGP %d %d %d %d %d %d %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]),
                            &(tmp_rgp[3]),
                             &(tmp_rgp[4]),
                             &(tmp_rgp[5]),
                             &(tmp_rgp[6]),
                             &(tmp_rgp[7]),
                             &(tmp_rgp[8]));
      if(records_read != 9) error_exit("Error parsing RGP number 4");
      if(tmp_rgp[0] != cc->num_rgroup) error_exit("Error RGroup number in num_rgp=4 \n");
      cc->name_rgroup[0] = tmp_rgp[2];  // Rgroup name is int format
      cc->name_rgroup[1] = tmp_rgp[4];
      cc->name_rgroup[2] = tmp_rgp[6];
      cc->name_rgroup[3] = tmp_rgp[8];

      cc->atomid_rgroup[0] = tmp_rgp[1];   //local atom index in mol
      cc->atomid_rgroup[1] = tmp_rgp[3];
      cc->atomid_rgroup[2] = tmp_rgp[5];
      cc->atomid_rgroup[3] = tmp_rgp[7];

      printf("RGP is %d %d %d %d %d %d %d %d %d\n", cc->num_rgroup,
             cc->atomid_rgroup[0], cc->name_rgroup[0],
             cc->atomid_rgroup[1], cc->name_rgroup[1],
             cc->atomid_rgroup[2], cc->name_rgroup[2],
             cc->atomid_rgroup[3], cc->name_rgroup[3]);
    }
    else{
      printf("number RGP is %d  allowed number is %d not supported\n", cc->num_rgroup, MAX_CONNECTIONS);

    }


    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF"); // M  END
    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF"); // $$$$

    records_read = sscanf(buffer, "%s", end_str);
    if(end_str[0] == '$'){
      printf("parser sdf finished\n");
    }
    fclose(fp);
  }
}

void read_sdf_v2000_withR(Sdf_MetaData *meta, BNG_Model *model) {
  int i;
  FILE *fp =NULL;
  char buffer[LINE_SIZE + 1];
  char *ptr;
  char end_str[4];
  int tmp;
  int records_read=0;

  Sdf_Ctab *cc=NULL;
  Sdf_Ctab_Counter *counter=NULL;

  int num_moltypes = model->num_moltypes;

  meta->num_data = num_moltypes;

  meta->ctabs = malloc(num_moltypes*sizeof(Sdf_Ctab));
  meta->counts = malloc(num_moltypes*sizeof(Sdf_Ctab_Counter));

  for(i=0; i < num_moltypes; i++){
    char filename[100];
    concatenateStrings(model->moltypes[i].name, "sdf", filename);
    fp = fopen(filename, "r");
    if(!fp){
      printf("***ERROR: error in %s file\n", filename);
      exit(EXIT_FAILURE);
    }

    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");

    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");

    cc = meta->ctabs+i;
    counter=meta->counts+i;

    records_read =sscanf(buffer, "%d %d %d %d %d %d %d %d %d %d %d %s\n",
                         &(counter->num_atoms), &(counter->num_bonds),
               &(counter->num_atom_lists), &tmp, /* for obsolete record */
               &(counter->flag_chiral), &(counter->num_stext_entries),
               &(counter->num_react_components), &(counter->num_reactants),
               &(counter->num_products), &(counter->num_intermediates),
               &(counter->num_additional_lines), counter->version);

    if (records_read != 10)
      error_exit("Error parsing sdf headline\n");

    counter->last_globalid = -1; // init
    cc->counter.num_atoms = counter->num_atoms;
    cc->counter.num_bonds = counter->num_bonds;

    int num_atoms = counter->num_atoms;
    int num_bonds = counter->num_bonds;
    printf(" index molec %d num_atoms is %d num_bonds is %d\n", i, num_atoms, num_bonds);

    cc->atoms = (Sdf_Ctab_Atom *)malloc(num_atoms*sizeof(Sdf_Ctab_Atom));
    cc->bonds = (Sdf_Ctab_Bond *)malloc(num_bonds*sizeof(Sdf_Ctab_Bond));

    int index_rg=0;
    Sdf_Ctab_Atom *atom=NULL;

    for (int n = 0; n < num_atoms; n++) {
      if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");
      atom = cc->atoms + n;
      ptr = buffer;
      while (isspace(*ptr)) ptr++;
      if (isdigit(*ptr)) {
        records_read =
            sscanf(ptr, "%lf %lf %lf %s %d %d %d %d %d %d %d %d %d %d %d %d",
                   &(atom->pos[0]), &(atom->pos[1]), &(atom->pos[2]), atom->symbol,
                   &(atom->mass_diff), &(atom->charge), &(atom->stereo_parity),
                   &(atom->hydrogen_count), &(atom->flag_stereo_care),
                   &(atom->valence), &(atom->h0_designator),
                   &(atom->type_react_component), &(atom->num_react_component),
                   &(atom->mapping_number), &(atom->flag_invertion),
                   &(atom->flag_exact_change));
        atom->global_id = -1;
        atom->atom_num = atsym2no(atom->symbol); //
        atom->mol_id = i; // biaoji

        if (records_read != 16) error_exit("Error parsing sdf ctab atom data");
      }
    }

    Sdf_Ctab_Bond *bond=NULL;

    for(int n=0; n< num_bonds; n++){
      bond = cc->bonds+n;
      if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");
      ptr = buffer;

      while(isspace(*ptr)) ptr++;
      if(isdigit(*ptr)){
        records_read = sscanf(ptr, "%d %d %d %d %d %d",
                              &(bond->fatom),
                              &(bond->satom),
                              &(bond->type),
                              &(bond->stereo),
                              &(bond->topology),
                              &(bond->reacting_center_status)
                              );
        if(records_read != 6) error_exit("Error parsing bond data");
      }
    }

    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");

    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF"); // M  END
    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF"); // $$$$

    records_read = sscanf(buffer, "%s", end_str);
    if(end_str[0] == '$'){
      printf("parser sdf finished\n");
    }
    fclose(fp);
  }
}


int read_pdb(char *filename, Sdf_Ctab *ctab){
  int num_atoms=0;
  double x,y,z;
  double tmpf1, tmpf2;
  int ielem=1;
  int i;
  int len;
  FILE *fp;
  char buffer[200];
  char *ptr;

  char tmpchar[6], tmpchar2[6];
  int atomindex[1];
  char atomname[6], atomname2[6];
  int tmpint;
  char elem[2];

  num_atoms = get_file_length(filename);
  fp = fopen(filename, "r");
  int records_read;
  if(!fp){
    printf("***ERROR: error in %s file\n", filename);
  }

  ctab->counter.num_atoms = num_atoms;
  ctab->atoms = (Sdf_Ctab_Atom *)malloc(sizeof(Sdf_Ctab_Atom)*num_atoms);

  Sdf_Ctab_Atom *atom=NULL;
  for(int i=0; i<num_atoms; i++){
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    atom = ctab->atoms + i;
    len=sscanf(buffer, "%s %d %s %s %s %d %lf %lf %lf %lf %lf\n",
               tmpchar, atomindex, atomname, atomname2, tmpchar2, &tmpint,
               &(atom->pos[0]), &(atom->pos[1]), &(atom->pos[2]), &tmpf1, &tmpf2);
    if(len !=11){
      error_exit("Error parsing atom hedl\n");
    }
    printf("tmpchar %s\n", tmpchar);
    printf("atomname %s atomname2 %s\n", atomname, atomname2);
    for(int j=0; j< strlen(atomname); j++){
      if(atomname[j]>='A' && atomname[j]<='Z' || atomname[j] >='a' && atomname[j]<='z')
        elem[j]= atomname[j];
      else if(atomname[j]>=0 && atomname[j] <=9){
        elem[j]='\0';
      }
    }
    atom->atom_num = atsym2no(elem);
    atom->global_id = -1;
    printf(" ctab->atoms %d %d %lf %lf %lf\n", i, atom->atom_num, atom->pos[0], atom->pos[1], atom->pos[2]);
  }
  fclose(fp);

}

/*
int main(int argc, const char **argv){

  BNG_Model *model=(BNG_Model*)malloc(sizeof(BNG_Model));
  alloc_BNG_Model(model, MAX_MOLECULES, MAX_SPECIES, MAX_REACTIONS, MAX_REACTIONS);
  read_bngl("in.bngl", model);
  print_model(model);

  Sdf_MetaData *meta=NULL;
  meta = (Sdf_MetaData *)malloc(sizeof(Sdf_MetaData));

  read_sdf_v2000_single(meta, model);

  free_BNG_Model(model);
  return 0;
}
*/
