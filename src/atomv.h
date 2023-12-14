/***********************************************************************
  atomv.h: an include file for atom.c
***********************************************************************/
#define Ratom_1 1.0         /* RGB color of an atom */
#define Gatom_1 0.2
#define Batom_1 0.2
#define Ratom_2 0.2
#define Gatom_2 0.2
#define Batom_2 1.0
#define MAX_LINE_LENGTH 100

typedef struct {          /* Atom data type */
  float crd[3];
} AtomType;

int nlon=18, nlat=9;      /* Number of polygons for a sphere in the 
                             longitudinal & lateral directions */
float atom_radius = 0.9;  /* Atomic radius in Lennard-Jones unit */
int winx=1280, winy=1280;   /* Window size */
float min_ext[3], max_ext[3];  
                          /* Range of atomic coordinates:
                             (left,lower,back), (right,top,front) */
int natoms;               /* number of atoms */
AtomType *atoms;          /* array of atoms */
int *atomType;
float eye[3];             /* position of eye point */
float center[3];          /* position of look reference point */
float up[3];              /* up direction for camera */

/**********************************************************************
Parameters to loop through the appropriate simulation files to animate
**********************************************************************/
char *basePat = "/Users/parkerh/Documents/research/USC/MXenes/dump/run_Ti3C2_450_270_40/dump.Ti3C2_450_270_4.0.";
int stepLim = 9000;
int stepSize = 1000;
int step = 0;

