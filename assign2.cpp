/*
  CSCI 480
  Assignment 2
 */

#include <stdio.h>
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <math.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "pic.h"



int g_iMenuId;

int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;

double s = 1.0/2;
//double basis[4][4];
double basisLater[16] = {-s, 2-s, s-2, s, s*2, s-3, 3-2*s, -s, -s, 0, s, 0, 0, 1, 0, 0};

GLuint texture[2];

double drawLine = 500.0;

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;

CONTROLSTATE g_ControlState = ROTATE;

/* state of the world */
float g_vLandRotate[3] = {0.0, 0.0, 0.0};
float g_vLandTranslate[3] = {0.0, 0.0, 0.0};
float g_vLandScale[3] = {1.0, 1.0, 1.0};


int oneThousand = 1000;
int totalPoint = 0;
int Uspeed = 0;

struct point * Position;
struct point * T;
struct point * B;
struct point * N;
struct point * TSpeedup;               //lists that hold the points, tangent, normal and binormal

double maxHeight = 0;

int printScreen = 0;
char myfile[2048];




/* represents one control point along the spline */
struct point {
   double x;
   double y;
   double z;
};

/* spline struct which contains how many control points, and an array of control points */
struct spline {
   int numControlPoints;
   struct point *points;
};

/* the spline array */
struct spline *g_Splines;

/* total number of splines */
int g_iNumOfSplines;

//print screen
void saveScreenshot (char *filename)
{
int i, j;   

Pic *in = NULL;   

Pic *out = NULL;

if (filename == NULL)       

return;

in = pic_alloc(640, 480, 3, NULL);   

out = pic_alloc(640, 480, 3, NULL);

printf("File to save to: %s\n", filename);
//  for (i=479; i>=0; i--) {

//    glReadPixels(0, 479-i, 640, 1, GL_RGB, GL_UNSIGNED_BYTE, &in->pix[i*in->nx*in->bpp]);

//  }

glReadPixels(0, 0, 640, 480, GL_RGB, GL_UNSIGNED_BYTE, &in->pix[0]);       

for ( int j=0; j<480; j++ ) { 

for ( int i=0; i<640; i++ ) { 

PIC_PIXEL(out, i, j, 0) = PIC_PIXEL(in, i, 480-1-j, 0); 

PIC_PIXEL(out, i, j, 1) = PIC_PIXEL(in, i, 480-1-j, 1);             

PIC_PIXEL(out, i, j, 2) = PIC_PIXEL(in, i, 480-1-j, 2); 

} 

}

if (jpeg_write(filename, out))       

printf("File saved Successfully\n");   

else       

printf("Error in Saving\n");

pic_free(in);    

pic_free(out);
}


// load the texture
void textLoad(int i, char *filename){
  Pic *img;
  img = jpeg_read(filename, NULL);
  glBindTexture(GL_TEXTURE_2D, texture[i]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img->nx, img->ny, 0, GL_RGB, GL_UNSIGNED_BYTE, &img->pix[0]);
  pic_free(img);
}

//matrix cross product
void crossProduct(struct point * a, struct point * b, struct point * result){
  result->x = a->y * b->z - a->z * b->y;
  result->y = a->z * b->x - a->x * b->z;
  result->z = a->x * b->y - a->y * b->x;

}

//matrix multiplicaiton
void matrixMultiply(double *a, double *b, double *result, int i, int j, int k){
  double s = 0;
      for(int c = 0; c < i; c++){
        for(int d = 0; d< k; d ++){
          for(int e = 0; e < j; e++){
            s = s + a[c * j + e] * b[e * k + d];
          }
          result[c * k + d] = s;
          s = 0;
        }
      }

}

// construct control points
void constructControlPoint(int i, double *control){
  struct point * points = g_Splines->points;
  for (int j = 0; j < 4; j++){

    control[j*3] = points[i+j].x;
    control[j*3 + 1] = points[i+j].y;
    control[j*3 + 2] = points[i+j].z;
  }
}

//calculate points
void pointCal(double u, double *control, struct point *result){
  double tempMatrix[4];
  tempMatrix[0] = u * u * u;
  tempMatrix[1] = u * u;
  tempMatrix[2] = u;
  tempMatrix[3] = 1;

  double temp1[4];
  double final[3];
  matrixMultiply(tempMatrix, basisLater,temp1, 1, 4, 4);

  matrixMultiply(temp1, control, final, 1, 4, 3);

  result->x = final[0];
  result->y = final[1];
  result->z = final[2];
}

//calculate tangent
void tangentCal(double u, double *control, struct point *result){
  double tempMatrix[4];
  tempMatrix[0] = 3 * u * u;
  tempMatrix[1] = 2 * u;
  tempMatrix[2] = 1;
  tempMatrix[3] = 0;

  double temp1[4];
  double final[3];
  matrixMultiply(tempMatrix, basisLater,temp1, 1, 4, 4);
  matrixMultiply(temp1, control, final, 1, 4, 3);
  result->x = final[0];
  result->y = final[1];
  result->z = final[2];
}

//turn to unit vector
void turnUnit(struct point * unit){

  double distance = sqrt(unit->x * unit->x + unit->y * unit->y + unit->z * unit->z);
  unit->x = unit->x / distance;
  unit->y = unit->y / distance;
  unit->z = unit->z / distance;
}

int loadSplines(char *argv) {
  char *cName = (char *)malloc(128 * sizeof(char));
  FILE *fileList;
  FILE *fileSpline;
  int iType, i = 0, j, iLength;


  /* load the track file */
  fileList = fopen(argv, "r");
  if (fileList == NULL) {
    printf ("can't open file\n");
    exit(1);
  }
  
  /* stores the number of splines in a global variable */
  fscanf(fileList, "%d", &g_iNumOfSplines);

  g_Splines = (struct spline *)malloc(g_iNumOfSplines * sizeof(struct spline));

  /* reads through the spline files */
  for (j = 0; j < g_iNumOfSplines; j++) {
    i = 0;
    fscanf(fileList, "%s", cName);
    fileSpline = fopen(cName, "r");

    if (fileSpline == NULL) {
      printf ("can't open file\n");
      exit(1);
    }

    /* gets length for spline file */
    fscanf(fileSpline, "%d %d", &iLength, &iType);

    /* allocate memory for all the points */
    g_Splines[j].points = (struct point *)malloc(iLength * sizeof(struct point));
    g_Splines[j].numControlPoints = iLength;

    /* saves the data to the struct */
    while (fscanf(fileSpline, "%lf %lf %lf", 
     &g_Splines[j].points[i].x, 
     &g_Splines[j].points[i].y, 
     &g_Splines[j].points[i].z) != EOF) {
      i++;
    }
  }

  free(cName);

  return 0;
}



//draw tracks
void drawSpline(){

    glColor3f(1.0, 1.0, 1.0);

//left rails
glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x - 0.03* B[i].x;
  temp.y = Position[i].y - 0.03* B[i].y;
  temp.z = Position[i].z - 0.03* B[i].z;
  temp1.x = Position[i+1].x - 0.03* B[i+1].x;
  temp1.y = Position[i+1].y - 0.03* B[i+1].y;
  temp1.z = Position[i+1].z - 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x + space * B[i].x + space * N[i].x, temp.y+ space * B[i].y + space * N[i].y, temp.z+ space * B[i].z + space * N[i].z);
  glVertex3f(temp.x + space * B[i].x - space * N[i].x, temp.y+ space * B[i].y - space * N[i].y, temp.z+ space * B[i].z - space * N[i].z);
  glVertex3f(temp1.x + space * B[i+1].x - space * N[i+1].x, temp1.y+ space * B[i+1].y - space * N[i+1].y, temp1.z+ space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp1.x + space * B[i+1].x + space * N[i+1].x, temp1.y+ space * B[i+1].y + space * N[i+1].y, temp1.z+ space * B[i+1].z + space * N[i+1].z);
  
}
glEnd();

glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x - 0.03* B[i].x;
  temp.y = Position[i].y - 0.03* B[i].y;
  temp.z = Position[i].z - 0.03* B[i].z;
  temp1.x = Position[i+1].x - 0.03* B[i+1].x;
  temp1.y = Position[i+1].y - 0.03* B[i+1].y;
  temp1.z = Position[i+1].z - 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x + space * B[i].x + space * N[i].x, temp.y+ space * B[i].y + space * N[i].y, temp.z+ space * B[i].z + space * N[i].z);
  glVertex3f(temp1.x + space * B[i+1].x + space * N[i+1].x, temp1.y+ space * B[i+1].y + space * N[i+1].y, temp1.z+ space * B[i+1].z + space * N[i+1].z);
  glVertex3f(temp1.x - space * B[i+1].x + space * N[i+1].x, temp1.y- space * B[i+1].y + space * N[i+1].y, temp1.z- space * B[i+1].z + space * N[i+1].z);
  glVertex3f(temp.x - space * B[i].x + space * N[i].x, temp.y- space * B[i].y + space * N[i].y, temp.z- space * B[i].z + space * N[i].z); 
  
}
glEnd();

glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x - 0.03* B[i].x;
  temp.y = Position[i].y - 0.03* B[i].y;
  temp.z = Position[i].z - 0.03* B[i].z;
  temp1.x = Position[i+1].x - 0.03* B[i+1].x;
  temp1.y = Position[i+1].y - 0.03* B[i+1].y;
  temp1.z = Position[i+1].z - 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x - space * B[i].x + space * N[i].x, temp.y- space * B[i].y + space * N[i].y, temp.z- space * B[i].z + space * N[i].z);
  glVertex3f(temp.x - space * B[i].x - space * N[i].x, temp.y- space * B[i].y - space * N[i].y, temp.z- space * B[i].z - space * N[i].z);
  glVertex3f(temp1.x - space * B[i+1].x - space * N[i+1].x, temp1.y- space * B[i+1].y - space * N[i+1].y, temp1.z- space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp1.x - space * B[i+1].x + space * N[i+1].x, temp1.y- space * B[i+1].y + space * N[i+1].y, temp1.z- space * B[i+1].z + space * N[i+1].z);
  
}
glEnd();

glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x - 0.03* B[i].x;
  temp.y = Position[i].y - 0.03* B[i].y;
  temp.z = Position[i].z - 0.03* B[i].z;
  temp1.x = Position[i+1].x - 0.03* B[i+1].x;
  temp1.y = Position[i+1].y - 0.03* B[i+1].y;
  temp1.z = Position[i+1].z - 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x + space * B[i].x - space * N[i].x, temp.y+ space * B[i].y - space * N[i].y, temp.z+ space * B[i].z - space * N[i].z);
  glVertex3f(temp1.x + space * B[i+1].x - space * N[i+1].x, temp1.y+ space * B[i+1].y - space * N[i+1].y, temp1.z+ space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp1.x - space * B[i+1].x - space * N[i+1].x, temp1.y- space * B[i+1].y - space * N[i+1].y, temp1.z- space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp.x - space * B[i].x- space * N[i].x, temp.y- space * B[i].y - space * N[i].y, temp.z- space * B[i].z - space * N[i].z); 
  
}
glEnd();





//right rails
glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x + 0.03* B[i].x;
  temp.y = Position[i].y + 0.03* B[i].y;
  temp.z = Position[i].z + 0.03* B[i].z;
  temp1.x = Position[i+1].x + 0.03* B[i+1].x;
  temp1.y = Position[i+1].y + 0.03* B[i+1].y;
  temp1.z = Position[i+1].z + 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x + space * B[i].x + space * N[i].x, temp.y+ space * B[i].y + space * N[i].y, temp.z+ space * B[i].z + space * N[i].z);
  glVertex3f(temp.x + space * B[i].x - space * N[i].x, temp.y+ space * B[i].y - space * N[i].y, temp.z+ space * B[i].z - space * N[i].z);
  glVertex3f(temp1.x + space * B[i+1].x - space * N[i+1].x, temp1.y+ space * B[i+1].y - space * N[i+1].y, temp1.z+ space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp1.x + space * B[i+1].x + space * N[i+1].x, temp1.y+ space * B[i+1].y + space * N[i+1].y, temp1.z+ space * B[i+1].z + space * N[i+1].z);
  
}
glEnd();

glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x + 0.03* B[i].x;
  temp.y = Position[i].y + 0.03* B[i].y;
  temp.z = Position[i].z + 0.03* B[i].z;
  temp1.x = Position[i+1].x + 0.03* B[i+1].x;
  temp1.y = Position[i+1].y + 0.03* B[i+1].y;
  temp1.z = Position[i+1].z + 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x + space * B[i].x + space * N[i].x, temp.y+ space * B[i].y + space * N[i].y, temp.z+ space * B[i].z + space * N[i].z);
  glVertex3f(temp1.x + space * B[i+1].x + space * N[i+1].x, temp1.y+ space * B[i+1].y + space * N[i+1].y, temp1.z+ space * B[i+1].z + space * N[i+1].z);
  glVertex3f(temp1.x - space * B[i+1].x + space * N[i+1].x, temp1.y- space * B[i+1].y + space * N[i+1].y, temp1.z- space * B[i+1].z + space * N[i+1].z);
  glVertex3f(temp.x - space * B[i].x + space * N[i].x, temp.y- space * B[i].y + space * N[i].y, temp.z- space * B[i].z + space * N[i].z); 
  
}
glEnd();

glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x + 0.03* B[i].x;
  temp.y = Position[i].y + 0.03* B[i].y;
  temp.z = Position[i].z + 0.03* B[i].z;
  temp1.x = Position[i+1].x + 0.03* B[i+1].x;
  temp1.y = Position[i+1].y + 0.03* B[i+1].y;
  temp1.z = Position[i+1].z + 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x - space * B[i].x + space * N[i].x, temp.y- space * B[i].y + space * N[i].y, temp.z- space * B[i].z + space * N[i].z);
  glVertex3f(temp.x - space * B[i].x - space * N[i].x, temp.y- space * B[i].y - space * N[i].y, temp.z- space * B[i].z - space * N[i].z);
  glVertex3f(temp1.x - space * B[i+1].x - space * N[i+1].x, temp1.y- space * B[i+1].y - space * N[i+1].y, temp1.z- space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp1.x - space * B[i+1].x + space * N[i+1].x, temp1.y- space * B[i+1].y + space * N[i+1].y, temp1.z- space * B[i+1].z + space * N[i+1].z);
  
}
glEnd();

glBegin(GL_QUADS);
for (int i = 0; i < totalPoint -1; i++){
  struct point  temp, temp1;
  temp.x = Position[i].x + 0.03* B[i].x;
  temp.y = Position[i].y + 0.03* B[i].y;
  temp.z = Position[i].z + 0.03* B[i].z;
  temp1.x = Position[i+1].x + 0.03* B[i+1].x;
  temp1.y = Position[i+1].y + 0.03* B[i+1].y;
  temp1.z = Position[i+1].z + 0.03* B[i+1].z;

  float space = 0.005;
  glVertex3f(temp.x + space * B[i].x - space * N[i].x, temp.y+ space * B[i].y - space * N[i].y, temp.z+ space * B[i].z - space * N[i].z);
  glVertex3f(temp1.x + space * B[i+1].x - space * N[i+1].x, temp1.y+ space * B[i+1].y - space * N[i+1].y, temp1.z+ space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp1.x - space * B[i+1].x - space * N[i+1].x, temp1.y- space * B[i+1].y - space * N[i+1].y, temp1.z- space * B[i+1].z - space * N[i+1].z);
  glVertex3f(temp.x - space * B[i].x- space * N[i].x, temp.y- space * B[i].y - space * N[i].y, temp.z- space * B[i].z - space * N[i].z); 
  
}
glEnd();


//rail cross-section
glColor3f(0.5f, 0.35f, 0.05f);
glBegin(GL_QUADS);
for (int i = 0; i < totalPoint ; i = i + 80){
    struct point  temp, temp1, temp2, temp3;
    int wide = 8;
  temp.x = Position[i].x + 0.03* B[i].x;
  temp.y = Position[i].y + 0.03* B[i].y;
  temp.z = Position[i].z + 0.03* B[i].z;
  temp1.x = Position[i+wide].x + 0.03* B[i+wide].x;
  temp1.y = Position[i+wide].y + 0.03* B[i+wide].y;
  temp1.z = Position[i+wide].z + 0.03* B[i+wide].z;
  temp2.x = Position[i].x - 0.03* B[i].x;
  temp2.y = Position[i].y - 0.03* B[i].y;
  temp2.z = Position[i].z - 0.03* B[i].z;
  temp3.x = Position[i+wide].x - 0.03* B[i+wide].x;
  temp3.y = Position[i+wide].y - 0.03* B[i+wide].y;
  temp3.z = Position[i+wide].z - 0.03* B[i+wide].z;

  glVertex3f(temp.x, temp.y, temp.z );
  glVertex3f(temp1.x, temp1.y, temp1.z );
  glVertex3f(temp3.x, temp3.y, temp3.z );
  glVertex3f(temp2.x, temp2.y, temp2.z );
}
glEnd();



}

//load textures
void loadScene(){

glBindTexture(GL_TEXTURE_2D, texture[0]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glEnable(GL_TEXTURE_2D);
  glBegin(GL_POLYGON);
  glColor3f(1.0, 1.0, 1.0); 
  glTexCoord2f(1.0, 0.0);
  glVertex3f(drawLine, drawLine, drawLine);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(-drawLine, drawLine, drawLine);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(-drawLine, drawLine, -drawLine);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(drawLine, drawLine, -drawLine);
  glEnd();

 glBindTexture(GL_TEXTURE_2D, texture[1]);
  glBegin(GL_POLYGON);
  glColor3f(1.0, 1.0, 1.0); 
  glTexCoord2f(1.0, 0.0);
  glVertex3f(drawLine, -drawLine, -drawLine);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(-drawLine, -drawLine, -drawLine);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(-drawLine, drawLine, -drawLine);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(drawLine,drawLine, -drawLine);
  glEnd();

 glBindTexture(GL_TEXTURE_2D, texture[1]);
  glBegin(GL_POLYGON);
  glColor3f(1.0, 1.0, 1.0); 
  glTexCoord2f(1.0, 0.0);
  glVertex3f(-drawLine, drawLine, drawLine);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(drawLine, drawLine, drawLine);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(drawLine, -drawLine, drawLine);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(-drawLine,-drawLine, drawLine);
  glEnd();

 glBindTexture(GL_TEXTURE_2D, texture[1]);
  glBegin(GL_POLYGON);
  glColor3f(1.0, 1.0, 1.0); 
  glTexCoord2f(1.0, 0.0);
  glVertex3f(-drawLine, drawLine, drawLine);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(-drawLine, -drawLine, drawLine);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(-drawLine, -drawLine, -drawLine);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(-drawLine,drawLine, -drawLine);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, texture[1]);
  glBegin(GL_POLYGON);
  glColor3f(1.0, 1.0, 1.0); 
  glTexCoord2f(1.0, 0.0);
  glVertex3f(drawLine, -drawLine, drawLine);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(drawLine, -drawLine, -drawLine);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(-drawLine, -drawLine, -drawLine);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(-drawLine,-drawLine, drawLine);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, texture[1]);
  glBegin(GL_POLYGON);
  glColor3f(1.0, 1.0, 1.0); 
  glTexCoord2f(1.0, 0.0);
  glVertex3f(drawLine, -drawLine, drawLine);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(drawLine, drawLine, drawLine);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(drawLine, drawLine, -drawLine);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(drawLine,-drawLine, -drawLine);
  glEnd();

  glDisable(GL_TEXTURE_2D);


}


  void myinit()
{
  /* setup gl view here */
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
 

  gluPerspective(60, 640.0/480, 0.01, 1000.0);
  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_MODELVIEW);
  glGenTextures(2, texture);
  textLoad(0, "ground.jpg");
  textLoad(1, "sky.jpg");


  /* calculate points' position, tangent, normal and binormal*/
  totalPoint = (g_Splines -> numControlPoints - 3) * (oneThousand+1);
  Position = (struct point*)malloc(totalPoint * sizeof(struct point));
  T = (struct point*)malloc(totalPoint * sizeof(struct point));
  TSpeedup = (struct point*)malloc(totalPoint * sizeof(struct point));
  B = (struct point*)malloc(totalPoint * sizeof(struct point));
  N = (struct point*)malloc(totalPoint * sizeof(struct point));
  struct point v = {0.0, 0.0, 1.0};
  double control[12];
  constructControlPoint(0, control);
  tangentCal(0, control, &T[0]);
  tangentCal(0, control, &TSpeedup[0]);
  turnUnit(&T[0]);
  
  crossProduct(&T[0], &v, &N[0]);
  turnUnit(&N[0]);
  crossProduct(&T[0], &N[0], &B[0]);
  turnUnit(&B[0]);

  for(int i = 0; i < g_Splines->numControlPoints - 3; i++){
    double control[12];
    constructControlPoint(i, control);
    for(int u = 0; u <=1000; u ++){
      pointCal(u*1.0/1000, control, &Position[1001 * i + u]);
      if(Position[1001 * i + u].z > maxHeight)
        maxHeight = Position[1001 * i + u].z;
      tangentCal(u*1.0/1000, control, &T[1001 * i + u]);
      tangentCal(u*1.0/1000, control, &TSpeedup[1001 * i + u]);
      turnUnit(&T[1001 * i + u]);
      if(i+ u > 0){
      crossProduct(&B[1001 * i + u - 1], &T[1001 * i + u], &N[1001 * i + u]);
      turnUnit(&N[1001 * i + u]);
      crossProduct(&T[1001 * i + u], &N[1001 * i + u], &B[1001 * i + u]);
      turnUnit(&B[1001 * i + u]);
      }
    }
  }


}



void display(){

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  glTranslatef(g_vLandTranslate[0], g_vLandTranslate[1], g_vLandTranslate[2]); //Translation

  glRotatef(g_vLandRotate[0], 1.0, 0.0, 0.0); 
  glRotatef(g_vLandRotate[1], 0.0, 1.0, 0.0);        // Rotation
  glRotatef(g_vLandRotate[2], 0.0, 0.0, 1.0);

  glScalef(g_vLandScale[0], g_vLandScale[1], g_vLandScale[2]);   // Scale

  /* set coordinates for the glulookat funcitons */
  struct point a,b;
  a.x = Position[Uspeed].x + 0.05 * N[Uspeed].x;
  a.y = Position[Uspeed].y + 0.05 * N[Uspeed].y;
  a.z = Position[Uspeed].z + 0.05 * N[Uspeed].z;

  b.x = Position[Uspeed].x + 0.05 * N[Uspeed].x + T[Uspeed].x;
  b.y = Position[Uspeed].y + 0.05 * N[Uspeed].y + T[Uspeed].y;
  b.z = Position[Uspeed].z + 0.05 * N[Uspeed].z + T[Uspeed].z;

  gluLookAt(a.x, a.y, a.z, b.x, b.y, b.z, N[Uspeed].x, N[Uspeed].y, N[Uspeed].z);


  loadScene();  

 drawSpline();

//set the speed close to real
float l = sqrt(pow(TSpeedup[Uspeed].x,2) + pow(TSpeedup[Uspeed].y,2) + pow(TSpeedup[Uspeed].z,2)) ;

  Uspeed = Uspeed+ 5 * sqrt(2* 9.8 * (maxHeight+ 0.1 - a.z))/l ;



  if(Uspeed >= totalPoint){
    Uspeed = 0;
  }


 
  glFlush();
  glutSwapBuffers();

}

void doIdle(){

  static int count = 0;
  if(printScreen == 1){
      if (count < 1000) {
        sprintf(myfile, "anim.%04d.jpg", count + 1);
        saveScreenshot(myfile);
        count++;
    }
// myFilenm will be anim.0001.jpg, anim.0002.jpg..........anim.0999.jpg // ..
}

  glutPostRedisplay();
}

void menufunc(int value)
{
  switch (value)
  {
    case 0:
      exit(0);
      break;
  }
}

void keyboard(unsigned char key, int x, int y){

  if(key =='d'||key =='D')
    printScreen = 1 - printScreen;  
}


void mousedrag(int x, int y)
{

  int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};
  
  switch (g_ControlState)
  {
    case TRANSLATE:  
      if (g_iLeftMouseButton)
      {
        g_vLandTranslate[0] += vMouseDelta[0]*0.01;
        g_vLandTranslate[1] -= vMouseDelta[1]*0.01;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandTranslate[2] += vMouseDelta[1]*0.01;
      }
      break;
    case ROTATE:
      if (g_iLeftMouseButton)
      {
        g_vLandRotate[0] += vMouseDelta[1];
        g_vLandRotate[1] += vMouseDelta[0];
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandRotate[2] += vMouseDelta[1];
      }
      break;
    case SCALE:
      if (g_iLeftMouseButton)
      {
        g_vLandScale[0] *= 1.0+vMouseDelta[0]*0.01;
        g_vLandScale[1] *= 1.0-vMouseDelta[1]*0.01;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandScale[2] *= 1.0-vMouseDelta[1]*0.01;
      }
      break;
  }
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mouseidle(int x, int y)
{
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mousebutton(int button, int state, int x, int y)
{

  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
 
  switch(glutGetModifiers())
  {
    case GLUT_ACTIVE_CTRL:
      g_ControlState = TRANSLATE;
      break;
    case GLUT_ACTIVE_SHIFT:
      g_ControlState = SCALE;
      break;
    default:
      g_ControlState = ROTATE;
      break;
  }

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}




int main (int argc, char ** argv)
{
  if (argc<2)
  {  
  printf ("usage: %s <trackfile>\n", argv[0]);
  exit(0);
  }
  loadSplines(argv[1]);
  glutInit(&argc,argv);
  


    glutInitDisplayMode (GLUT_DOUBLE |GLUT_DEPTH |GLUT_RGBA);

    glutInitWindowSize(640,480);
    glutInitWindowPosition(0,0);
    glutCreateWindow("Assign2");
    //glutReshapeFunc();

    glEnable(GL_DEPTH_TEST); // enable depth test

  glutDisplayFunc(display);

  /* allow the user to quit using the right mouse button menu */
  g_iMenuId = glutCreateMenu(menufunc);
  glutSetMenu(g_iMenuId);
  glutAddMenuEntry("Quit",0);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mousedrag);
  /* callback for idle mouse movement */
  glutPassiveMotionFunc(mouseidle);
  /* callback for mouse button changes */
  glutMouseFunc(mousebutton);

  glutKeyboardFunc(keyboard);
  


  myinit();
  glutMainLoop();

  return 0;
}
