// https://github.com/StructuralBioinformaticsLab/libmol/blob/6d1c3f158a7dccadf6ee60098f1f1dbaf423be22/mol.0.0.6/sdf.c
#ifndef DEFINE_H
#define DEFINE_H

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#define LINELEN_MAX 160

// AmberTools_version.update_number
//const char *const ANTECHAMBER_VERSION="24.0" ;

#define COLORTEXT "YES"
#define REDUCEFACTOR 1
#define MAXCHAR 2048
#define MAXATOM 25600
#define MAXRES  5120
#define MAXBOND 51200
#define MAXRING 5000
#define MAXGAS  2000            /*maximum gasiteger parameters */
/*
For MAXRING, no dynamic memory applied since the actual number is determined
using a recursive function. However, for small and middle-sized molecules,
it is unlikely that the ring num is larger than 1000
*/
#define MAXCYCLE 1000
#define OUTPUTSTEP 10
#define MAXTWIST 10
#define ECSLONG 2
#define COSCUT 120
#define DEGRAD 3.1415926/180
#define VDWIDIST 10
#define ESIDIST 14
#define THETACUT 15
#define CUBE 2.0
#define MAXWILDATOM 20
#define MAXSCHAIN 100
#define MAXCES 20
#define MAXBEED 20
#define MAXATOMTYPE 250
#ifdef NCSU_PENALTIES_H
# define PSCUTOFF 10
#else
# define PSCUTOFF 10
# define MAXVASTATE 8192
# endif                         /* NCSU_PENALTIES_H */
#define MAX_CES_BOND 100
#define MAXAPUNIT 20
#define MAXAPCONS 20
#define MAXBLFPARM 1024


#define RAD_TO_DEG 57.2957795131
#define DEG_TO_RAD 0.0174532925199
#define degreesToRadians(angleDegrees) (angleDegrees * PI / 180.0)
#define radiansToDegrees(angleRadians) (angleRadians * 180.0 / PI)
#define EPSILON 8.8817841970012523e-016 /* 4.0 * DBL_EPSILON */


#define PI 3.141592653589793
#define TWOPI 2 * PI
// try max times to gen position
#define GENPOS_TRY_MAX 50


// temp
#define P_MAX_EL 118

static int periodic_max_el=P_MAX_EL;


static char *table[P_MAX_EL+2]={"X","H","He",
  "Li","Be", "B", "C", "N", "O", "F","Ne",
  "Na","Mg","Al","Si", "P", "S","Cl","Ar",
  "K","Ca",
  "Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
  "Ga","Ge","As","Se","Br","Kr",
  "Rb","Sr",
   "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
  "In","Sn","Sb","Te", "I","Xe",
  "Cs","Ba",
  "La",
  "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
  "Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn",
  "Fr","Ra",
  "Ac",
  "Th","Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lw",
  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
0};


/*
* corresponding table of VDW radii. g03
*/

static  double espfitvdwr[102] = {0.0, 1.2, 1.2,                   //None H He
                           1.37, 1.45, 1.45, 1.5, 1.5, // Li Be B C N
                           1.4, 1.35, 1.3,             // O F Ne
                           1.57, 1.36, 1.24, 1.17, 1.8,// Na Mg Al Si P
                           1.75, 1.7, 1.88,            // S Cl Ar
                           2.75, 0.0,                  // K Ca(na)
                          0.0, 0.0, 0.0, 0.0, 1.44,     // Sc Ti V Cr Mn
                          1.43, 0.0, 1.63, 1.4, 1.39,   // Fe Co Ni Cu Zn
                          1.87, 1.86, 2.0, 2.0, 2.3,  // Ga Ge As Se Br
                          2.02,                        // Kr
                          0.0, 0.0,                    // Rb Sr
                          0.0, 0.0, 0.0, 0.0, 0.0,     // Y Zr Nb Mo Tc
                          0.0, 0.0, 1.63, 1.72, 1.58,  // Ru Rh Pd Ag Cd
                          1.93, 2.17, 2.2, 2.2, 2.15,  // In Sn Sb Te I
                          2.16,                        // Xe
                          0.0, 0.0,                    // Cs Ba
                          0.0, 0.0, 0.0, 0.0, 0.0,     // La Ce Pr Nd Pm
                          0.0, 0.0, 0.0, 0.0, 0.0,     // Sm Eu Gd Tb Dy
                          0.0, 0.0, 0.0, 0.0, 0.0,      // Ho Er Tm Yb Lu
                          0.0,                         // Hf
                          0.0, 0.0, 0.0, 0.0, 0.0,     // Ta W Re Os Ir
                          1.72, 1.66, 1.55, 1.96, 1.02,// Pt Au Hg Tl Pb
                          0.0, 0.0, 0.0, 0.0, 0.0,     // Bi Po At Rn Fr
                          0.0, 0.0, 0.0,               // Ra Ac Th
                          0.0, 1.86, 0.0, 0.0, 0.0,     // Pa U Np Pu Am
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0};                       // Md

static double covalent_radius[102] = {
    0.1,  0.32, 0.93,                                // Non H He
    1.23, 0.90, 0.82, 0.77, 0.75, 0.73, 0.72, 0.71,  // Li Be B C N 3  O F Ne
    1.54, 1.36, 1.18, 1.11, 1.06, 1.02, 0.99, 0.98,  // Na Mg Al Si P S Cl Ar
    2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17,  // K Ca Sr Ti V Cr Mn Fe
    1.16, 1.15, 1.17, 1.25, 1.26, 1.22, 1.20, 1.16,  // Co Ni Cu Zn Ga Ge As Se
    1.14, 1.12, 2.16, 1.91, 1.62, 1.45, 1.34, 1.30,  // Br Kr Rb Sr Y Zr Nb Mo
    1.27, 1.25, 1.25, 1.28, 1.34, 1.48, 1.44, 1.41,  // Tc Ru Rh Pd Ag Cd In Sn
    1.40, 1.36, 1.33, 1.31, 2.35, 1.98, 1.69, 1.65,  // Sb Te I Xe Cs Ba La Ce
    1.65, 1.64, 1.63, 1.62, 1.85, 1.61, 1.59, 1.59,  // Pr Nd Pm Sm Eu Gd Tb Dy
    1.58, 1.57, 1.56, 1.74, 1.56, 1.44, 1.34, 1.30,  // Ho Er Tm Yb Lu Hf Ta W
    1.28, 1.26, 1.27, 1.30, 1.34, 1.49, 1.48, 1.47,  // Re Os Ir Pt Au Hg Tl Pb
    1.46, 1.46, 1.45, 1.50, 2.60, 2.21, 2.15, 2.06,
    2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 0.1, 0.1, 0.1, 0.1,
};

static const double an2masses[119] = {
0.0,1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406,
12.0,14.00307400478,15.99491461956,18.998403224,19.99244017542,
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629,
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983,
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292,
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676,
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932,
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297,
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757,
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089,
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109,
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011,
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271,
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325,
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108,
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724,
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500,
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451,
285.183698,287.191186,292.199786,291.206564,293.214670};


#define NEIGHBOR_MAX 20

// define max bond number of each node
#define MAX_CONNECTIONS 4
#define MAX_BONDS 100

#define LINE_SIZE 300

#define LOG_ERROR(f, ...) \
  fprintf(stderr, "%s:%d: " f "\n", __FILE__, __LINE__, __VA_ARGS__)


#define SD_SDF_DELIMITER "$$$$"
#define SDF_ATOMS_NUM_MAX 255
#define SDF_BONDS_NUM_MAX 255
#define SDF_CHIRAL 1
#define SDF_NOT_CHIRAL 0
#define SDF_CURRENT_CTAB_VER "V2000"
#define SDF_ATOM_LIST "L"
#define SDF_ATOM_UNSPECIFIED "*"
#define SDF_ATOM_R_GROUP "R#"

#define SDF_NON_STEREO_PARITY 0
#define SDF_EVEN_STEREO_PARITY 1
#define SDF_ODD_STEREO_PARITY 2
#define SDF_EITHER_OR_UNMARKED_STEREO_CENTER 3

#define SDF_IGNORE_STEREO_CONFIG 0
#define SDF_STEREO_CONFIG_MUST_MATCH 1
#define SDF_REACTANT 1
#define SDF_PRODUCT 2
#define SDF_INTERMEDIATE 3
#define SDF_CONFIGURATION_INVERTED 1
#define SDF_CONFIGURATION_RETAINED 2

#define SDF_SINGLE_BOND 1
#define SDF_DOUBLE_BOND 2
#define SDF_TRIPLE_BOND 3
#define SDF_AROMATIC_BOND 4
#define SDF_SINGLE_OR_DOUBLE_BOND 5
#define SDF_SINGLE_OR_AROMATIC_BOND 6
#define SDF_DOUBLE_OR_AROMATIC_BOND 7
#define SDF_ANY_BOND 8

#define SDF_RING_BOND 1
#define SDF_CHAIN_BOND 2

#define SDF_CTAB_COUNTS_RECORDS_NUM 12

// ?
#define TOL  0.1
#define ARR3SET(a,b) {*a = *b; *(a+1) = *(b+1); *(a+2) = *(b+2);}


typedef struct{
  char name[10];
  int id;
} BNG_Component;

typedef struct {
  char name[50];
  int num_comp;
  BNG_Component *components;
} BNG_MoleculeType;

/*
  not use now
typedef struct{
  char name[50];
  char lable[50];
  int *edges;
  int *adjacency;
}BNG_SpeciesGraph;
*/
typedef struct {
  char name[50];
  int count;
  //  BNG_SpeciesGraph *graph;
} BNG_Species;

typedef struct{
  char reactants[100];
  char products[100];
  double rate;
} BNG_ReactionRule;

typedef struct{
  char aname[100];
  char bname[100];
  int bid[2];
  double rate;
} BNG_BindingRule;

typedef struct{
  int num_moltypes;
  int num_species;
  int num_rules;
  int num_bindings;
  BNG_MoleculeType *moltypes;
  BNG_Species *species;
  BNG_ReactionRule *rules;
  BNG_BindingRule *bindings;
} BNG_Model;


#define MAX_LINE_LENGTH 100
#define MAX_MOLECULES 100
#define MAX_SPECIES 100
#define MAX_REACTIONS 100
#define MAX_COMPONENTS 5

typedef struct{
  int num_mols;
  int max_attempts;
  double cutoff_scalar;
  double site_cutoff;
  double pos_cutoff;
  char sdf_file[100];
  char cyc_file[100];
  int *cyc;
  int num_cycle;
 } Settings;

//void read_settings(Settings *set);
//void init_settings(Settings *set);
//void free_settings(Settings *set);
//void print_settings(Settings *set);

typedef struct{
  int num_atoms;
  int num_bonds;
  int last_globalid;
  int num_atom_lists; // max 30?
  int flag_chiral;
  int num_stext_entries;
  int num_react_components;
  int num_reactants;
  int num_products;
  int num_intermediates;
  int num_additional_lines;
  char version[16]; // V2000
} Sdf_Ctab_Counter;

typedef struct{
  char symbol[4];
  double pos[3];
  double charge;
  int mass_diff;
  int stereo_parity;
  int hydrogen_count;
  int flag_stereo_care;
  int valence;
  int type_react_component;
  int mapping_number;
  int flag_invertion;
  int flag_exact_change;
  int h0_designator;
  int num_react_component;
  int global_id;
  int atom_num;
  int mol_id;
  char atomtype[4];
  char resname[4];
} Sdf_Ctab_Atom;

typedef struct{
  int fatom; // first_atom_number
  int satom; // second_atom_number
  int type;
  int stereo;
  int topology;
  int reacting_center_status;
  Sdf_Ctab_Atom sca;
  Sdf_Ctab_Atom scb;
} Sdf_Ctab_Bond;


typedef struct {
  Sdf_Ctab_Counter counter;
  Sdf_Ctab_Atom *atoms;
  Sdf_Ctab_Bond *bonds;
  int num_rgroup;
  int *name_rgroup; // save each rgoup name id //
  int *atomid_rgroup; // rgroup's site in mol
} Sdf_Ctab;

typedef struct {
  Sdf_Ctab_Counter *counts; // counts[num_data]
  Sdf_Ctab *ctabs; // num_data ctabs[num_data]
  int num_data;
} Sdf_MetaData;

typedef struct{
  int num_mol;
  int total_num;
  int **groupid;   // save rgroup site [nmol][max_bond] == groupid
  int *rgp; // save each rgroup name id
  int *num_rgp; // save each mol's num rg
  int *name;
  int **atomid;
} RgroupRecord;

// simple define
typedef struct {
  int num_atoms;
  int *names; // connect to Sdf_Ctab_Atom's global_id
  double *coords;
  int *atom_nums;
} Molecule;

// geometry_utils
// 3 point present a geometry plane
typedef struct{
  double a1[3];
  double a2[3];
  double a3[3];
}Plane;

/* A pointer to filehandle and its real name */
/* Used for user defined file IO operations */
struct file_stream {
  char *name;   /* File name */
  FILE *stream; /* File handle structure */
};
/* Symbol hash table entry */
/* Used to parse and store user defined symbols from the MDL input file */
struct sym_entry {
  struct sym_entry *next; /* Chain to next symbol in this bin of the hash */
  int sym_type;           /* Symbol Type */
  char *name;             /* Name of symbol*/
  void *value;            /* Stored value, cast by sym_type */
  int count;
};

/* Symbol hash table */
/* Used to parse and store user defined symbols from the MDL input file */
struct sym_table_head {
  struct sym_entry **entries;
  int n_entries;
  int n_bins;
};
/* Linked list of symbols */
/* Used to parse and retrieve user defined symbols having wildcards from the
 * MDL input file */
struct sym_table_list {
  struct sym_table_list *next;
  struct sym_entry *node; /* Symbol table entry matching a user input wildcard
                             string */
};

/* Doubly linked list of object names */
struct name_list {
  struct name_list *next;
  struct name_list *prev;
  char *name; /* An object name */
};



typedef struct {
    char name[10];
    char aa[20];
    char ambername[10];
    char bcctype[10];
    char element[10];
    char resid[10];
    char chain[5];
    double x;
    double y;
    double z;
    double radius;
    double charge;
    double eps;
    double lja;
    double ljb;
    double bond;
    double angle;
    double twist;
    int id;
    int ter;
    int resno;
    int resseqno;
    int hdonor;
    int haccept;
    int atomtype;
    int con[6];
    int bondatom;
    int angleatom;
    int twistatom;
    int connum;
    int mainatom;
    int atomicnum;
    int select;
    int type;
    int type2;                  /* 1: rotatable bond atoms, 0: irrotatable bond atoms */
    int type3;                  /* 1: ring atoms, 0: non-ring atoms */
    int ewd;                    /* electrwithdraw or not: 1 ewd; -1 edonor; 0 others */
    int arom;                   /* 5: 5 mem.ring aromatic atom; 6: 6 mem.ring aromatic atom;

                                   -1: aliphatic; 0: others */
    int aliph;
    int saturate;               /* 1: saturated; -1 unsaturated; */
    int improper;
    int ifreeze;    /* freeze the torsional angle when this atom is the fourth atom of the torsional angle*/
} ATOM;
typedef struct {
    char name[10];
    double x;
    double y;
    double z;
} ATOMSML;
typedef struct {
    char name[10];
    int atomtype;
    int con[6];
    double x;
    double y;
    double z;
    double charge;
} ATOMDBS;
typedef struct {
    char name[10];
    double x;
    double y;
    double z;
    double radius;
    double charge;
    double eps;
    double lja;
    double ljb;
    int hdonor;
    int haccept;
    int atomtype;
    int con[6];
} LIGATOM;
typedef struct {
    int atm1;
    int atm2;
    int atm3;
    int atm4;
    char type;
    double k;
    double kexi;
    int multi;
} TWIST;
typedef struct {
    double x;
    double y;
    double z;
} POINT;
typedef struct {
    int atomno;
    double x;
    double y;
    double z;
    double esp;
    double dist;
} ESP;
typedef struct {
    double x;
    double y;
    double z;
    double total;
} DM;
typedef struct {
    double xx;
    double yy;
    double zz;
    double xy;
    double xz;
    double yz;
} QM;
typedef struct {
    int atconnum;
    int apindex;
    char ap[512];
    char atname[10];
    char cesname[10];
} CHEMENV;

typedef struct {
    int bondi;
    int bondj;
    int type;
    int type2;                  /* 0- irrotatable bond, 1- rotatable bond */
    int jflag;                  /* 1- already assigned bond type, 0 or -1 no */
    double bcc;
} BOND;

typedef struct {
    char name1[10];
    char name2[10];
    double length;
    double force;
    int iprint;
    int attn;
} BOND_FF;

typedef struct {
    int num;
    int atomno[12];
} RING;

typedef struct {
    int nr;
    int ar1;
    int ar2;
    int ar3;
    int ar4;
    int ar5;
    int rg[12];
} AROM;

typedef struct {
    char dkeyword[MAXCHAR];     /*divcon keyword */
    char ekeyword[MAXCHAR];     /*empirical calculation keyword, copied from dkeyword or mkeyword */
    char gkeyword[MAXCHAR];     /*gaussian keyword */
    char gm[MAXCHAR];           /*gaussian %mem */
    char gdsk[MAXCHAR];         /*gaussian %disk */
    char gn[MAXCHAR];           /*gaussian %nproc */
    char gesp[MAXCHAR];         /*gaussian esp file */
    char mkeyword[MAXCHAR];     /*mopac keyword */
    char skeyword[MAXCHAR];     /*sqm keyword */
    char gopt[MAXCHAR]; /*gaussian keyword for optimization*/
    char gsp[MAXCHAR];  /*gaussian keyword for single point calculation*/
    char tor[MAXCHAR];  /*torsinal angles*/
    char resname[MAXCHAR];
    char chkfile[20];
    char atom_type_def[10];
    char resfilename[MAXCHAR];
    char gfilename[MAXCHAR];
    char connect_file[MAXCHAR];
    char radius_file[MAXCHAR];
    char longresname[MAXCHAR];
    char inf_filename[MAXCHAR];
    int divcon;                 /* 1-use divcon, 0-not use divcon */
    int multiplicity;           /*molecular multiplicity */
    int icharge;                /*integer charge */
    int usercharge;             /*user read in charge */
    int igkeyword;              /*user provides keyword or not */
    int igopt;  /*user gopt keyword or not*/
    int igsp;  /*user gsp keyword or not*/
    int itor;  /*user tor keyword or not*/
    int igdsk;  /*user gdsk keyword or not*/
    double dcharge;             /*float charge */
    int gv;                     /*gaussian version flag: 1- g09, 0- other versions */
    int eqcharge;
} MOLINFO;

typedef struct {
    char intype[10];
    char outtype[10];
    char atype[10];
    char chargetype[50];
    int rnindex;                /*read in residue index ? */
    int intstatus;              /*information status */
    int pfindex;                /*purify intermediate files */
    int prediction_index;       /*the index of performing atomtype */
    int bpindex;                /*bond type prediciton index */
    int maxatom;
    int maxbond;
    int maxring;
    int max_path_length;        /*maximum path length to judge atom equilibration */
    int verify_pdb_atomname;    /*varify pdb atom names using the possible element field in a pdb file */
    int  atseq;     /*atomic sequence changeable? 1 for yes or 0 for no */
} CONTROLINFO;

typedef struct {
    char name[10];
    double a;
    double b;
    double c;
    double d;
    double charge;
} GASTEIGER;

typedef struct {
    char elem[10];
    int atomicnum;
    int flag;                   /* = 0, use the default radii in Gaussian */
    int bs;                     /* basis set ID */
    double vdw;
    double mk;
} ESPPARM;

typedef struct {
    int id;
    char bs[MAXCHAR];
} BASISSET;

typedef struct {
    char name[10];
} NAME;

typedef struct {
    int aps[10];
} AV;

/*parameters for calculating bond stretching force constant*/
typedef struct {
    double bfkai;
    double refbondlength;
    int id1;
    int id2;
    char elem1[10];
    char elem2[10];
} BLF_PARM;

/*parameters for calculating bond bending force constant*/
typedef struct {
    double anglec;
    double anglez;
    int id;
    char elem[10];
} BAF_PARM;
//BLF_PARM blf_parm[MAXBLFPARM];
//BAF_PARM baf_parm[120];
//double blf_exp_const = 4.5;     /*exponential parameter for calculating blf */




static inline int rotl(const int x, int k) {
  return (x << k) | (x >> (64 - k));
}

static inline uint64_t  splitmix64( uint64_t x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}


/// Sets the interal state of the random number generator.
static inline void rng_seed( uint64_t *s, uint64_t seed){
  s[0] = splitmix64(seed);
  s[1] = splitmix64(s[0]);
}

/// Gives a random long integer between 0 and 2^(64)-1
static inline uint64_t rng_next( uint64_t *s) {
  const uint64_t s0 = s[0];
  uint64_t s1 = s[1];
  const uint64_t result = s0 + s1;

  s1 ^= s0;
  s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
  s[1] = rotl(s1, 37); // c

  return result;
}


/// Gives a random double uniformly distributed in the
/// range (0,1]:  exclusive zero, inclusive one
static inline double rng_uniform( uint64_t *s){
  uint64_t rand_int = rng_next(s);
  uint64_t bits = 0x3ff0000000000000 | (rand_int >> 12);
  return 2.0 - *((double *) &bits);
}



static char* atno2sym(unsigned no){
  if (no>P_MAX_EL) return (0);
  else return (table[no]);
}

static unsigned int atsym2no(char* sym){
  char *p1,*p2;
  int i;

  for(i=1;i<P_MAX_EL;i++){
    p1=sym;
    p2=table[i];
    while(*p1==' ') p1++;
    while((*p1)&&(*p2)&&(!(((*p1)^(*p2))&95))){p1++;p2++;}
    if(((*p1==0)||(*p1==' ')||(*p1==':'))&&(*p2==0)) return i;
  }

  return(0);
}


static double get_pte_vdw(const int idx)
{
  if((idx < 1) || (idx >=112)) return espfitvdwr[0];
  return espfitvdwr[idx];
}
// get covlent radius
static double get_pte_cov(const int idx)
{
  if((idx < 1) || (idx >=112)) return espfitvdwr[0];
  return covalent_radius[idx];
}


static void error_exit(char *msg){
  fprintf(stderr, "Aborting: %s\n", msg);
  exit(1);
}


static inline void vector3_subtract(double a[3], double b[3], double c[3]){
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
  return;
}


static inline void vector3_mat3b3_multiply(double a[3][3], double b[3],
                                           double c[3]) {
  double temp[3];
  temp[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
  temp[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
  temp[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
  c[0] = temp[0];
  c[1] = temp[1];
  c[2] = temp[2];
  return;
}


static inline void vector3_add(double a[3], double b[3], double sum[3]) {
  sum[0] = a[0] + b[0];
  sum[1] = a[1] + b[1];
  sum[2] = a[2] + b[2];
  return;
}

#define ARR3MUL(a,b) {*a *= *b; *(a+1) *= *(b+1); *(a+2) *= *(b+2);}
#define VEC3MUL(v,m){*v *= m; *(v+1) *= m; *(v+2) *= m;}

static inline void print_settings(Settings *set){
  printf("INPUT SETTINGS:\n");
  printf("-----------------------------\n");
  printf("Input Sdf2000 file path is: %s\n", set->sdf_file);
  printf("Number of molecule to link(same with sdf_file):         %d \n", set->num_mols);
  printf("Maximum attempts per molecule is:       %d\n", set->max_attempts);
  //  printf("Random seed:  %d\n", set->random_seed);
  printf("molecule cutoff scalar is %f\n", set->cutoff_scalar);
  printf("site cutoff is %f\n", set->site_cutoff);
  printf("next pos distance cutoff is %f\n", set->pos_cutoff);
}

static inline void init_settings(Settings *set){
  set->cutoff_scalar =0.65;
  set->pos_cutoff = 1.6;
  set->site_cutoff = 1.3;
  set->max_attempts =500;
  set->num_cycle=0;
  set->cyc = NULL;
}


static inline void free_settings(Settings *set) {
  free(set->cyc);
  free(set);
}

static inline void read_settings(Settings *set){
  FILE *fileptr;
  size_t len = 0;
  char *line = NULL;
  char *sub_line = NULL;
  int read;
  char sdf_file[100];
  char cyc_file[100];
  fileptr = fopen("control.in", "r");
  if (!fileptr) {
    printf("***ERROR: no control.in file \n");
    exit(EXIT_FAILURE);
  }

  // read from control
  while ((read = getline(&line, &len, fileptr)) != -1) {

    // if comment
    if (strstr(line, "#") != NULL)
      continue;

    sub_line = strtok(line, " ");

    if (strcmp(sub_line, "sdf_file") == 0) {
      sub_line = strtok(NULL, " ");
      strcpy(sdf_file, sub_line);
      sdf_file[strcspn(sdf_file, "\n")] = 0;
      strcpy(set->sdf_file, sdf_file);
      continue;
    }
    if (strcmp(sub_line, "cyc_file") == 0) {
      sub_line = strtok(NULL, " ");
      strcpy(cyc_file, sub_line);
      cyc_file[strcspn(cyc_file, "\n")] = 0;
      strcpy(set->cyc_file, cyc_file);
      continue;
    }

    if (strcmp(sub_line, "num_mols") == 0) {
        sub_line = strtok(NULL, " ");
        set->num_mols = atoi(sub_line);
        continue;
      }

    if (strcmp(sub_line, "cutoff_scalar") == 0) {
        sub_line = strtok(NULL, " ");
        set->cutoff_scalar = atof(sub_line);
        continue;
      }

    if (strcmp(sub_line, "pos_distance_cutoff") == 0) {
        sub_line = strtok(NULL, " ");
        set->pos_cutoff = atof(sub_line);
        continue;
      }

    if (strcmp(sub_line, "site_distance_cutoff") == 0) {
        sub_line = strtok(NULL, " ");
        set->site_cutoff = atof(sub_line);
        continue;
      }

    if (strcmp(sub_line, "max_attempts") == 0) {
        sub_line = strtok(NULL, " ");
        set->max_attempts= atoi(sub_line);
        continue;
      }
  }
  fclose(fileptr);

}

#endif
