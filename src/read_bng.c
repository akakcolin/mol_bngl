#define PCRE2_CODE_UNIT_WIDTH 8
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pcre2.h"
#include <stdbool.h>
#include "define.h"


void alloc_BNG_Model(BNG_Model *model, int num_mol, int num_species, int num_rules, int num_bindings){
  int i;
  model->num_moltypes = num_mol;
  model->num_species = num_species;
  model->num_rules = num_rules;
  model->num_bindings = num_bindings;

  model->moltypes=(BNG_MoleculeType *)malloc(sizeof(BNG_MoleculeType)*num_mol);

  for(i=0; i<model->num_moltypes;i++){
    model->moltypes[i].components = (BNG_Component *)malloc(sizeof(BNG_Component)*MAX_COMPONENTS);
    model->moltypes[i].num_comp = MAX_COMPONENTS;
  }
  model->species=(BNG_Species *)malloc(sizeof(BNG_Species)*num_species);
  model->rules=(BNG_ReactionRule *)malloc(sizeof(BNG_ReactionRule)*num_rules);
  model->bindings=(BNG_BindingRule *)malloc(sizeof(BNG_BindingRule)*num_bindings);

}

void free_BNG_Model(BNG_Model *model){
  int i;
  for(i=0; i<model->num_moltypes;i++){
    free(model->moltypes[i].components);
  }
  free(model->moltypes);
  free(model->species);
  free(model->rules);
  free(model->bindings);
  free(model);
}

void print_model(BNG_Model *model){
  for(int i=0; i<model->num_moltypes; i++){
    printf("%s \n", model->moltypes[i].name);
    for(int j=0; j<model->moltypes[i].num_comp; j++){
      printf(" (%s %d)\n", model->moltypes[i].components[j].name, model->moltypes[i].components[j].id);
    }
  }

  for(int i=0; i<model->num_species; i++){
    printf("%s %d\n", model->species[i].name, model->species[i].count);
  }
  for(int i=0; i<model->num_bindings; i++){
    printf("%s %s rate is %lf\n", model->bindings[i].aname, model->bindings[i].bname, model->bindings[i].rate);
  }

}

int parse_molecule_types(char *sources, BNG_MoleculeType *mol){
  pcre2_code *re_code, *re_code2, *re_code3;
  PCRE2_SPTR pattern, pattern2, pattern3;
  PCRE2_SPTR subject, subject2, subject3;
  int errnum, errnum2, errnum3;
  int rc, rc2, rc3;
  size_t subject_length, subject2_length, subject3_length;

  PCRE2_SIZE erroffset, erroffset2, erroffset3;
  PCRE2_SIZE *ovector, *ovector2, *ovector3;
  pcre2_match_data *match_data, *match_data2, *match_data3;
  char components[100];

  char *re_left="([A-Za-z0-9_]+)\\(([^)]+)\\)";
  char *re_comp="([^!]+)";
  char *re_split="([^, ]+)";

  pattern=(PCRE2_SPTR)re_left;
  pattern2=(PCRE2_SPTR)re_split;
  pattern3=(PCRE2_SPTR)re_comp;

  // use to molecule
  re_code = pcre2_compile(pattern, PCRE2_ZERO_TERMINATED, 0, &errnum, &erroffset, NULL);  // use default
  if (re_code == NULL) {
    PCRE2_UCHAR buffer[256];
    pcre2_get_error_message(errnum, buffer, sizeof(buffer));
    printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroffset,
           buffer);
    return 1;
  }
  match_data = pcre2_match_data_create_from_pattern(re_code, NULL);
  // use to component
  re_code2 =
      pcre2_compile(pattern2, PCRE2_ZERO_TERMINATED, 0, &errnum2, &erroffset2,
                    NULL);  // use default
  if (re_code2 == NULL) {
    PCRE2_UCHAR buffer[256];
    pcre2_get_error_message(errnum2, buffer, sizeof(buffer));
    printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroffset2,
           buffer);
    return 1;
  }
  match_data2 = pcre2_match_data_create_from_pattern(re_code2, NULL);

  // use to component
  re_code3 = pcre2_compile(pattern3,
                           PCRE2_ZERO_TERMINATED, 0, &errnum3, &erroffset3,
                           NULL);  // use default
  if (re_code3 == NULL) {
    PCRE2_UCHAR buffer[256];
    pcre2_get_error_message(errnum3, buffer, sizeof(buffer));
    printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroffset3,
           buffer);
    return 1;
  }
  match_data3 = pcre2_match_data_create_from_pattern(re_code3, NULL);

  // first subject
  subject = (PCRE2_SPTR)sources;
  subject_length = strlen((char *)subject);

  rc = pcre2_match(re_code, subject, subject_length, 0, 0, match_data, NULL);
  if (rc < 0) {
    switch (rc) {
      case PCRE2_ERROR_NOMATCH:
        printf("No Match\n");
        break;
      default:
        printf("Matching error %d\n", rc);
        break;
    }
    pcre2_match_data_free(match_data);
    pcre2_code_free(re_code);
    return 1;
  }
  ovector = pcre2_get_ovector_pointer(match_data);
  PCRE2_SPTR substring_start;
  size_t substring_length;
  if(rc==3){
    substring_start = subject + ovector[2];
    substring_length = ovector[3] - ovector[2];
    //printf("%2d: %.*s\n", 1, (int)substring_length, (char *)substring_start);
    strncpy(mol->name, (char *)substring_start, (int)substring_length);
    //printf("mol->name is %s\n", mol->name);

    substring_start = subject + ovector[4];
    substring_length = ovector[5] - ovector[4];
    //printf("%2d: %.*s\n", 2, (int)substring_length, (char *)substring_start);
    strncpy(components, (char *)substring_start, (int)substring_length);
  }

  pcre2_match_data_free(match_data);
  pcre2_code_free(re_code);


  subject2=(PCRE2_SPTR)components;
  subject2_length = strlen((char *)subject2);

  rc2 = pcre2_match(re_code2, subject2, subject2_length, 0, 0, match_data2, NULL);
  if (rc2 < 0) {
    switch (rc2) {
      case PCRE2_ERROR_NOMATCH:
        printf("No Match\n");
        break;
      default:
        printf("Matching error %d\n", rc2);
        break;
    }
    pcre2_match_data_free(match_data2);
    pcre2_code_free(re_code2);
    return 1;
  }
  ovector2 = pcre2_get_ovector_pointer(match_data2);

  int count=0;

  //mol->num_comp=rc2-1;

  BNG_Component *comp_buffer=(BNG_Component *)malloc(sizeof(BNG_Component)*MAX_COMPONENTS);

  while(rc2 >0){
    PCRE2_SIZE start_offset = pcre2_get_ovector_pointer(match_data2)[0];
    PCRE2_SIZE end_offset = pcre2_get_ovector_pointer(match_data2)[1];
    PCRE2_SPTR substring = subject2 + start_offset;
    size_t substring_length = end_offset - start_offset;
    //printf("component Match %d: %.*s\n", count + 1, (int)substring_length, substring);
    strncpy(comp_buffer[count].name, (char *)substring, (int)substring_length);
    rc2 = pcre2_match(re_code2, subject2, subject2_length, end_offset, 0, match_data2, NULL);
    count++;
  }
  mol->num_comp = count;

  mol->components=(BNG_Component *)malloc(sizeof(BNG_Component)*count);

  pcre2_match_data_free(match_data2);
  pcre2_code_free(re_code2);

  for (int i = 0; i < count; i++) {
    subject3 = (PCRE2_SPTR)comp_buffer[i].name;
    subject3_length = strlen((char *)subject3);
    //printf("%d %s\n", i, comp_buffer[i].name);
    rc3 = pcre2_match(re_code3, subject3, subject3_length, 0, 0, match_data3, NULL);
    if (rc3 < 0) {
      switch (rc3) {
        case PCRE2_ERROR_NOMATCH:
          printf("No Match\n");
          break;
        default:
          printf("Matching error %d\n", rc3);
          break;
      }
      pcre2_match_data_free(match_data3);
      pcre2_code_free(re_code3);
      return 1;
    }
    ovector3 = pcre2_get_ovector_pointer(match_data3);

    if (rc3 == 2) {
      PCRE2_SIZE start_offset = pcre2_get_ovector_pointer(match_data3)[0];
      PCRE2_SIZE end_offset = pcre2_get_ovector_pointer(match_data3)[1];
      PCRE2_SPTR substring = subject3 + start_offset;
      size_t substring_length = end_offset - start_offset;
      //printf("single Match %d: %.*s\n",  1, (int)substring_length, substring);
      strncpy(mol->components[i].name, (char *)substring, (int)substring_length);
      //      printf("%s\n", mol->components[i].name);
      rc3 = pcre2_match(re_code3, subject3, subject3_length, end_offset, 0,
                        match_data3, NULL);
      start_offset = pcre2_get_ovector_pointer(match_data3)[0];
      end_offset = pcre2_get_ovector_pointer(match_data3)[1];
      substring = subject3 + start_offset;
      mol->components[i].id = atoi(substring);
      //printf("%d\n", mol->components[i].id);
    }
  }
  free(comp_buffer);
  pcre2_match_data_free(match_data3);
  pcre2_code_free(re_code3);
  return 0;
}

int parse_species(char *sources, BNG_Species *specie){
  char *token = strtok(sources, " ");
  strcpy(specie->name, token);
  token = strtok(NULL, " ");
  specie->count = atoi(token);
}

int parse_binding(char *sources, BNG_BindingRule *rule){
  char *token = strtok(sources, " ");
  int field = 0;
  while (token != NULL && field < 5) {
    switch (field) {
      case 0:  // aname
        strncpy(rule->aname, token, sizeof(rule->aname) - 1);
        break;
      case 1:  // aid
        rule->bid[0] = atoi(token);
        break;
      case 2:  // bname
        strncpy(rule->bname, token, sizeof(rule->bname) - 1);
        break;
      case 3:  // bid
        rule->bid[1] = atoi(token);
        break;
      case 4:  // rate
        rule->rate = atof(token);
        break;
    }
    token = strtok(NULL, " ");
    field++;
  }
}

int read_bngl(char *filename, BNG_Model *model){
  FILE *fileptr;
  size_t len = 0;
  char *line = NULL;
  char *sub_line = NULL;
  int read;

  int mol_count=0;
  int species_count=0;
  int reaction_count=0;
  int in_mol_types=0;
  int in_species=0;
  int in_reaction_rules=0;
  int in_binding_rules=0;
  int in_param=0;
  int num_mol_types=0;
  int num_species=0;
  int num_rules=0;
  int num_bindings=0;
  int skip_line=1;
  int blank=0;

  fileptr = fopen(filename, "r");

  while((read = getline(&line, &len, fileptr)) != -1){
    if(strstr(line, "#") !=NULL)
      continue;
    if(read == 1 && line[0] == '\n')
      continue;

    sub_line = strtok(line, "\n");

    if(strcmp(sub_line, "begin molecule types") ==0){
      in_mol_types=1;
      skip_line=1;
      in_species =0;
      in_param =0;
      in_reaction_rules=0;
      printf("in_mol\n");
      continue;
    }else{
      skip_line=0;
    }

    if(strcmp(sub_line, "end molecule types") ==0){
      in_mol_types=0;
      continue;
    }

    if(strcmp(sub_line, "begin seed species") ==0){
      skip_line = 1;
      in_species=1;
      //      printf("in_speced\n");

      in_mol_types=0;
      in_param=0;
      in_reaction_rules =0;
      continue;
    }else{
      skip_line = 0;
    }

    if(strcmp(sub_line, "end seed species") ==0){
      in_species=0;
      continue;
    }

    if(strcmp(sub_line, "begin reaction rules") ==0){
      in_reaction_rules=1;
      skip_line = 1;
      continue;
    }else{
      skip_line=0;
    }

    if(strcmp(sub_line, "end reaction rules") ==0){
      in_reaction_rules=0;
      continue;
    }


    if(strcmp(sub_line, "begin binding") ==0){
      in_binding_rules=1;
      skip_line = 1;
      continue;
    }else{
      skip_line=0;
    }

    if(strcmp(sub_line, "end binding") ==0){
      in_binding_rules=0;
      continue;
    }


    if(in_mol_types == 1 && skip_line == 0){
      parse_molecule_types(sub_line, &model->moltypes[num_mol_types]);
      num_mol_types++;
    }
    if(in_species==1 && skip_line == 0){
      parse_species(sub_line, &model->species[num_species]);
      num_species++;
    }
    if(in_reaction_rules==1 && skip_line ==0){
      num_rules++;
    }

    if(in_binding_rules==1 && skip_line ==0){
      parse_binding(sub_line, &model->bindings[num_bindings]);
      num_bindings++;
    }
  }

  model->num_moltypes = num_mol_types;
  model->num_species = num_species;
  model->num_rules = num_rules;
  model->num_bindings = num_bindings;
  fclose(fileptr);

}


int main(int argc, char *argv[])
{
  BNG_Model *model=(BNG_Model*)malloc(sizeof(BNG_Model));
  alloc_BNG_Model(model, MAX_MOLECULES, MAX_SPECIES, MAX_REACTIONS, MAX_REACTIONS);
  read_bngl("in.bngl", model);
  print_model(model);
  //model->moltypes = (BNG_MoleculeType*)malloc(sizeof(BNG_MoleculeType));
  //model->moltypes->components=NULL;

  //char *sources="_ISOPS_(cooh!3, cooh!4)";

  // model->num_moltypes=1;

  print_model(model);

  free_BNG_Model(model);

  return 0;
}
