#ifndef _MORPHING_H
#define _MORPHING_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "libmesh5.h"
#include "memory.h"
#include "sparse.h"

#define LS_VER   "2.0a"
#define LS_REL   "Oct, 2017"
#define LS_CPY   "Copyright (c) ISCD "
#define LS_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

#define LS_Ver       (1<<0)
#define LS_Edg       (1<<1)
#define LS_Tri       (1<<2)
#define LS_Tet       (1<<3)

#define LS_LAMBDA     10.0e6
#define LS_MU         8.2e6
#define LS_E          1e14
#define LS_NU         0.05
#define LS_MAT        50
#define LS_CL         50
#define LS_RES        1.e-6
#define LS_MAXIT      5000
#define LS_TGV        1.e+30

#define LS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define LS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )

#define LONMAX   8192
#define PRECI    1.0
#define EPS      1.e-6
#define EPS1     1.e-20
#define EPS2     1.e-10
#define EPSD     1.e-30
#define EPSA     1.e-200
#define EPST    -1.e-2
#define EPSR     1.e+2
#define INIVAL_2d     1 //0.25
#define INIVAL_3d     3.0

#define MESH_DREF 2
#define MESH_ELREF 2
#define MESH_BREF 99

#define DEFAULT_NIT 400; // default iterations number
#define DEFAULT_TOL 1e-3; // default error allowed
#define DEFAULT_SAVE 10; // default intermediate step for saving


enum {None=0, Dirichlet, Load};
enum {P1=0, P2};

typedef struct {
    double   c[3];
    int      ref,s,mark;
} Point;
typedef Point * pPoint;

typedef struct {
    int     v[2],ref;
} Edge;
typedef Edge * pEdge;

typedef struct {
    int     v[3],ref,mark;
    double    g[3];
} Tria;
typedef Tria * pTria;

typedef struct {
    double    g[3];
    int     v[4],s,ref,mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
    double   delta,err,gr[3],tol;
    int      ncpu,nit,save_it;
    char     imprim,typ,cg,rhs;
   ;
} Info;

typedef struct {
    double  lambda,mu;
    int     ref;
} Mat;
typedef Mat * pMat;

typedef struct {
    double   u[3];
    int      ref;
    char     typ,elt,att;
} Cl;
typedef Cl * pCl;


typedef struct {
    int      np,na,nt,ne,hmax,hcur,ver,dim,mark,ref1[10],ref2[10],ref3[10],nref1,nref2,nref3,npbn,refb;
    int      *adja,*tab;
    char     *name,cltyp;
    double   min[3],max[3],delta;

    pCl      cl;
    pMat     mat;
    pPoint   point;
    pEdge    edge;
    pTria    tria;
    pTetra   tetra;
} Mesh;
typedef Mesh * pMesh;


typedef struct {
    int      dim,ver,np,nit,iter,size[2],nmat,nbcl;
    double  *u,*p,*d,err,*valp1;
    char    *name,cltyp;

    pMat     mat;
    pCl      cl;

} Sol;
typedef Sol * pSol;


typedef struct {
    int   ia,ib,k,nxt;
} hedge;

typedef struct {
    int     siz,max,nxt;
    hedge  *item;
} Hash;

typedef struct {
    int     size;
    int    *head;
    int    *link;
} Bucket;
typedef Bucket * pBucket;

typedef struct {
    Bucket bucket;
    Sol sol_distance;
    Sol sol_omega0;
    Mesh mesh_distance;
    Mesh mesh_omega0;
    Csr *A;
} Instance;
typedef Instance * pInstance;


/* prototypes */
int     loadMesh(pMesh );
int     loadSol(pSol );
int     saveSol(pSol, pMesh, pMesh, int );
int     saveContour(pMesh, pMesh);
int     scaleMesh(pMesh );
int     hashar_2d(pMesh ,Hash,int);
int     hashar_3d(pMesh,Hash,int);
int     hashP2(Hash ,int ,int );
pCl     getCl(pSol ,int,int);
int     getMat(pSol,int,double *,double *);
int     saveMesh(pMesh, pMesh, int);
int     moveMesh(pMesh,pSol,double);
int     invmatg(double[9],double[9]);
double  volume(double *,double *,double *,double *);
double  area_3d(double *,double *,double *,double *);
int     getRefel(pMesh, int);
int     getRefb(pMesh, int);
int     getRefm(pMesh, int);
int     elas(pInstance, double *);
int     minimization(pInstance);
int     saveDistance(pSol , pMesh , pMesh, int);
double  hausdorff_2d(pMesh, pMesh);
double  hausdorff_3d(pMesh, pMesh);
double  errdist_2d(pMesh , pMesh );
double  errdist_3d(pMesh , pMesh );
double  distpt_3d(pPoint p0,pPoint p1,pPoint p2,pPoint pq,char *proj);
double  distpt_23d(pPoint p1,pPoint p2,pPoint pa);
int     findstep(pInstance, double *);
int     computedistance(pInstance);
double  errL2_2d(pInstance);
double  errL2_3d(pInstance );

int     initMesh_3d(pMesh mesh, double c[3], double);
int     initMesh_2d(pMesh mesh, double c[2], double);



int     mshin1(pMesh mesh1,pSol sol1,pMesh mesh2,pSol sol2);
int     locateTetra(pMesh mesh,int nsdep,int base,double *p,double *cb);
int     locateTria(pMesh mesh,int nsdep,int base,double *p,double *cb);

pBucket newBucket_2d(pMesh ,int );
pBucket newBucket_3d(pMesh ,int );
int     buckin_2d(pMesh ,pBucket ,double *);
int     buckin_3d(pMesh ,pBucket ,double *);
int     intpp0_3d(pMesh ,pSol ,double *,int ,double *,double );
int     intpp0_2d(pMesh ,pSol ,double *,int ,double *,double );
int     intpp1_2d(pSol ,int *,double *,int ,double *);
int     intpp1_3d(pSol ,int *,double *,int ,double *);
int     locelt_2d(pMesh ,int ,double *,double *);
int     locelt_3d(pMesh ,int ,double *,double *);
int     closept_2d(pMesh ,double *);
int     closept_3d(pMesh ,double *);
int     hashelt_3d(pMesh );
int     hashelt_2d(pMesh );
int     boulep_2d(pMesh ,int ,int ,int *);
int     boulep_3d(pMesh ,int ,int ,int *);

/* function pointers */
pCsr    matA_P1_3d(pMesh,pSol);
pCsr    matA_P1_2d(pMesh,pSol);
pCsr    (*matA_P1)(pMesh,pSol);

int     evalderFunctional3D(pInstance, double *);
int     evalderFunctional2D(pInstance, double *);
int     (*evalderFunctional)(pInstance, double *);

double  evalFunctional3D(pInstance, double);
double  evalFunctional2D(pInstance, double);
double  (*evalFunctional)(pInstance, double);


pBucket (*newBucket)(pMesh ,int );
int     (*buckin)(pMesh ,pBucket ,double *);
int     (*locelt)(pMesh ,int ,double *,double *);
int     (*closept)(pMesh ,double *);
int     (*intpp0)(pMesh ,pSol ,double *,int ,double *,double);
int     (*intpp1)(pSol ,int *,double *,int ,double *);
int     (*hashelt)(pMesh );
int     (*boulep)(pMesh ,int ,int ,int *);


#endif
