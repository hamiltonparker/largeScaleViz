/***********************************************************************
  Program atomv.c--ball representation of atoms.
  Required files
    atomv.h:   Include file
    md.conf:   MD configuration file containing atomic coordinates
***********************************************************************/
#include "atomv.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>    /* Header file for the OpenGL library */
#include <OpenGL/glu.h>   /* Header file for the GLu library */
#include <GLUT/glut.h>    /* Header file for the GLut library */

GLuint sphereid;          /* display-list id of atom sphere geom */
GLuint atomsid;           /* display-list id of all atoms */
GLdouble fovy, aspect, near_clip, far_clip;  
                          /* parameters for gluPerspective() */
FILE *fp;                 /* pointer to open an MD-configuration file */

int nprocs, myid;

/* Function prototypes ************************************************/
void reshape(int, int);
void makeFastNiceSphere(GLuint, double);
void makeAtoms(void);
void makeCurframeGeom(void);
void drawScene(void);
void display(void);
void initView(float *, float *);
void readConf(void);
void animate(void);
void writeFrame(void);
void PPMWriter(unsigned char *, char*, int, int);
void setStepRange(int, int);

/**********************************************************************/
void reshape (int w, int h) {
/***********************************************************************
  Callback for glutReshapeFunc()
***********************************************************************/
  /* set the GL viewport to match the full size of the window */
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  aspect = w/(float)h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fovy,aspect,near_clip,far_clip);
  glMatrixMode(GL_MODELVIEW);
}

/**********************************************************************/
void makeFastNiceSphere(GLuint listid, double radius) {
/***********************************************************************
Called once to generate and compile sphere geometry into the given
display list id.
***********************************************************************/
  int i,j;
  float lon,lat;
  float loninc,latinc;
  float x,y,z;

  loninc = 2*M_PI/nlon;
  latinc = M_PI/nlat;

  glNewList(listid,GL_COMPILE);

    /* South-pole triangular fan */
    glBegin(GL_TRIANGLE_FAN);
      glNormal3f(0,-1,0);
      glVertex3f(0,-radius,0);
      lon = 0;
      lat = -M_PI/2 + latinc;
      y = sin(lat);
      for (i=0; i<=nlon; i++) {
        x = cos(lon)*cos(lat);
        z = -sin(lon)*cos(lat);
        glNormal3f(x,y,z);
        glVertex3f(x*radius,y*radius,z*radius);
        lon += loninc;
      }
    glEnd();

    /* Quadrilateral stripes to cover the sphere */
    for (j=1; j<nlat-1; j++) {
      lon = 0;
      glBegin(GL_QUAD_STRIP);
        for (i=0; i<=nlon; i++) {
          x = cos(lon)*cos(lat);
          y = sin(lat);
          z = -sin(lon)*cos(lat);
          glNormal3f(x,y,z);
          glVertex3f(x*radius,y*radius,z*radius);
          x = cos(lon)*cos(lat+latinc);
          y = sin(lat+latinc);
          z = -sin(lon)*cos(lat+latinc);
          glNormal3f(x,y,z);
          glVertex3f(x*radius,y*radius,z*radius);
          lon += loninc;
        }
      glEnd();
      lat += latinc;
    }

    /* North-pole triangular fan */
    glBegin(GL_TRIANGLE_FAN);
      glNormal3f(0,1,0);
      glVertex3f(0,radius,0);
      y = sin(lat);
      lon = 0;
      for (i=0; i<=nlon; i++) {
        x = cos(lon)*cos(lat);
        z = -sin(lon)*cos(lat);
        glNormal3f(x,y,z);
        glVertex3f(x*radius,y*radius,z*radius);
        lon += loninc;
      }
    glEnd();

  glEndList();
}

/**********************************************************************/
void makeAtoms() {
/***********************************************************************
  Makes display-list of all atoms in the current frame using spheres.
***********************************************************************/
  int i;
  float rval1,gval1,bval1;
  float rval2,gval2,bval2;

  glNewList(atomsid, GL_COMPILE);
  rval1 = Ratom_1; gval1 = Gatom_1; bval1 = Batom_1;  /* RGB color of an atom */
  rval2 = Ratom_2; gval2 = Gatom_2; bval2 = Batom_2;
  for (i=0; i < natoms; i++) {
    glPushMatrix();
    glTranslatef(atoms[i].crd[0],atoms[i].crd[1],atoms[i].crd[2]);
    if (atomType[i] == 1){
      glColor3f(rval1,gval1,bval1);
    }
    else if (atomType[i] == 2){
      glColor3f(rval2,gval2,bval2);
    }
    glCallList(sphereid);
    glPopMatrix();
  }
  glEndList();
}

/**********************************************************************/
void makeCurframeGeom() {
/***********************************************************************
  Reads the atoms information for the current time frame and makes the
  display-list of all the atoms' geometry.
***********************************************************************/
  makeAtoms();
}

/**********************************************************************/
void drawScene() {
/***********************************************************************
  Called by display() to draw the view of the current scene.
***********************************************************************/
  /* Define viewing transformation */
  gluLookAt(
    (GLdouble)eye[0],(GLdouble)eye[1],(GLdouble)eye[2],
    (GLdouble)center[0],(GLdouble)center[1],(GLdouble)center[2],
    (GLdouble)up[0],(GLdouble)up[1],(GLdouble)up[2]);
  glCallList(atomsid);
}

/**********************************************************************/
void display() {
/***********************************************************************
  Callback for glutDisplayFunc().  It clears the frame and depth 
  buffers and draws the atoms in the current frame.
***********************************************************************/
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  drawScene();
  glutSwapBuffers();
}

/**********************************************************************/
void initView (float *min_ext, float *max_ext) {
/***********************************************************************
  Initializes global viewing, lighting, and projection values.
***********************************************************************/
  GLfloat light_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_position1[] = {0.3, 0.5, 1.0, 0.0};
  GLfloat light_specular[]  = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_position2[] = {0.3, 0.5, 1.0, 0.0};
  float dif_ext[3],dis;
  int i;

  /* Define normal light */
  glLightfv(GL_LIGHT0,GL_DIFFUSE,light_diffuse);
  glLightfv(GL_LIGHT0,GL_POSITION,light_position1);
  glLightfv(GL_LIGHT1,GL_SPECULAR,light_specular);
  glLightfv(GL_LIGHT1,GL_POSITION,light_position2);

  /* Enable a single OpenGL light */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  /* Use depth buffering for hidden surface elimination */
  glEnable(GL_DEPTH_TEST);

  /* get diagonal and average distance of extent */
  for (i=0; i<3; i++) dif_ext[i] = max_ext[i] - min_ext[i];
  dis = 0.0;
  for (i=0; i<3; i++) dis += dif_ext[i]*dif_ext[i];
  dis = (float)sqrt((double)dis);

  /* set center in world space */
  for (i=0; i<3; i++) center[i] = min_ext[i] + dif_ext[i]/2.0;

  /* set initial eye & look at location in world space */
  eye[0] = center[0] - dis*0.25*0.6;
  eye[1] = center[1];
  eye[2] = center[2] + dis*0.4*0.6;
  up[0] = 1.0;
  up[1] = 0.0;
  up[2] = 1.0;

  /* set parameters for gluPerspective() */
  /* Near- & far clip-plane distances */
  near_clip = (GLdouble)( 0.1*(dis-0.5*dif_ext[2]) );
  far_clip  = (GLdouble)( 2.0*(dis+0.5*dif_ext[2]) );
  /* Field of view */
  fovy = (GLdouble)( 0.5*dif_ext[1]/(dis-0.5*dif_ext[2]) );
  fovy = (GLdouble)( 2*atan((double)fovy)/M_PI*180.0 );
  fovy = (GLdouble)(1.2*fovy);

  /* Enable the color material mode */
  glEnable(GL_COLOR_MATERIAL);
}

/**********************************************************************/
void readConf() {
/***********************************************************************
Read atomic coordinates from an MD-configuration file & allocates 
necessary arrays.
***********************************************************************/
  int l, j;

  /* Open an MD-configuration file */
  fp = fopen("md.conf","r");
  /* Read the # of atoms */
  fscanf(fp,"%d",&natoms);
  /* allocate atoms array */
  atoms = (AtomType *) malloc(sizeof(AtomType)*natoms);
  /* Maximum & minimum extent of system in angstroms */
  for (l=0; l<3; l++) fscanf(fp,"%f%f",&min_ext[l],&max_ext[l]);
  /* Atomic coordinates */
  for (j=0; j<natoms; j++)
    fscanf(fp,"%f %f %f",&(atoms[j].crd[0]),&(atoms[j].crd[1]),
                         &(atoms[j].crd[2]));
  fclose(fp);
}

void readDump(char* pattern) {
/*
Read atomic coordinates from lammps dump file
*/
  int l, j;
  int isAtoms = 0;
  char line[MAX_LINE_LENGTH] = {0};
  FILE * file = fopen(pattern,"r");

  while (fgets(line, MAX_LINE_LENGTH, file)){
      if (strstr(line, "ITEM: BOX BOUNDS") != NULL) {
        for (l=0; l<3; l++) fscanf(file,"%f%f",&min_ext[l],&max_ext[l]);
        continue;
      }
      if (strstr(line, "ITEM: NUMBER OF ATOMS") != NULL){
        fscanf(file,"%d",&natoms);
        continue;
      }
      if (strstr(line, "ITEM: ATOMS") != NULL){
        break;
      }
  }

  atoms = (AtomType *) malloc(sizeof(AtomType)*natoms);
  atomType = (int *) malloc(sizeof(int)*natoms);
  for (j=0; j<natoms; j++) {
    fgets(line, MAX_LINE_LENGTH, file);
    sscanf(line, "%*d %d %f %f %f %*f %*f %*f ",&(atomType[j]), &(atoms[j].crd[0]),
                                &(atoms[j].crd[1]),&(atoms[j].crd[2]));
  }
  fclose(fp);

}

void animate() {

  if (step < stepLim) {
    char snum[12] = {0};
    sprintf(snum, "%d", step);
    char fPat[strlen(basePat)+10];
    strcpy(fPat,basePat);
    strcat(fPat,snum);
    readDump(fPat);
    makeCurframeGeom();
    glutPostRedisplay();
    writeFrame();
    step = step + stepSize;
  }

}

void writeFrame() {
  unsigned char *image = (unsigned char*) malloc(sizeof(unsigned char)*3*winx*winy);
  glReadPixels(0, 0, winx, winy, GL_RGB, GL_UNSIGNED_BYTE, image);

  char buffer [33];
  sprintf(buffer,"../frames/img_frame_%d.ppm",step);
  PPMWriter(image,buffer,winx,winy);
}

void PPMWriter(unsigned char *in, char *name, int dimx, int dimy) {
  int i, j;
  FILE *file = fopen(name,"w");
  (void) fprintf(file, "P3 \n%d %d \n255\n",dimx, dimy);
  
  size_t tmp;

  for (j = 0; j < dimy; j++) {
    for (i = 0; i < dimx; i++) {
      static unsigned char color[3];
      color[0] = in[3*i+3*(dimy-j-1)*dimy];
      color[1] = in[3*i+3*(dimy-j-1)*dimy+1];
      color[2] = in[3*i+3*(dimy-j-1)*dimy+2];
      fprintf(file, "%d %d %d ", color[0],color[1],color[2]);
    }
    fprintf(file,"\n");
  }
  (void) fclose(file);
}

void setStepRange(int nprocs, int myid) {
  int stepRange = stepLim - step;
  int nSteps = stepRange/stepSize/nprocs;

  step = myid * nSteps;
  if (myid != nprocs-1){
    stepLim = nSteps*(myid + 1);
  }

}

/**********************************************************************/
int main(int argc, char **argv) {
/**********************************************************************/

  MPI_Init(&argc, &argv);
  
  glutInit(&argc, argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  setStepRange(nprocs,myid);

  /* Read atomic coordinates from an MD-configuration file */
  readDump("/Users/parkerh/Documents/research/USC/MXenes/dump/run_Ti3C2_450_270_40/dump.Ti3C2_450_270_6.0.0");

  /* Set up an window */
  /* Initialize display mode */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  /* Specify window size */
  glutInitWindowSize(winx, winy);
  /* Open window */
  glutCreateWindow("Lennard-Jones Atoms");

  /* Initialize view */
  initView(min_ext, max_ext);

  /* Set a glut callback functions */
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(animate);

  /* generate an OpenGL display list for single sphere */
  sphereid = glGenLists(1);
  makeFastNiceSphere(sphereid,atom_radius);
  
  /* generate an OpenGL display list for the atoms' geometry */
  atomsid = glGenLists(1);
  /* make the geometry of the current frame's atoms */
  makeCurframeGeom();

  /* Start main display loop */
  glutMainLoop();
  
  MPI_Finalize();
  return 0;
}
/**********************************************************************/
