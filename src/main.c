#include "morphing.h"
#define BUCKSIZ 64 //BUCKSIZ defines the discretization of cube/square [0,1[^dim
                   // in BUCKSIZ^dim points (Warning: mesh must be in [0,1[^dim)

Info info;

static void excfun(int sigid) {
    fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
    switch(sigid) {
        case SIGABRT:
            fprintf(stdout,"  Abnormal stop\n");  break;
        case SIGBUS:
            fprintf(stdout,"  Code error...\n");  break;
        case SIGFPE:
            fprintf(stdout,"  Floating-point exception\n"); break;
        case SIGILL:
            fprintf(stdout,"  Illegal instruction\n"); break;
        case SIGSEGV:
            fprintf(stdout,"  Segmentation fault.\n");  break;
        case SIGTERM:
        case SIGINT:
            fprintf(stdout,"  Programm killed.\n");  break;
    }
    fprintf(stdout," No data file saved.\n");
    exit(1);
}

static void usage(char *prog) {
    fprintf(stdout,"\n usage: %s [opts..] target[.mesh] template[.mesh]\n",prog);

    fprintf(stdout,"\n** Generic options :\n");
    fprintf(stdout,"-h       Print this message\n");

    fprintf(stdout,"\n** Parameters\n");
    fprintf(stdout,"-dref   nref [refs...]   reference of the Dirichlet boundary ( default 1 %d ) \n",MESH_DREF);
    fprintf(stdout,"-elref  nref [refs...]   reference of the Dirichlet elements ( default 1 %d ) \n",MESH_ELREF);
    fprintf(stdout,"-bref   nref [refs...]   reference of the surface to not be clamped or subjected to loads ( default 1 %d ) \n",MESH_BREF);
    fprintf(stdout,"-nit    n                iterations (default %d) \n",info.nit);
    fprintf(stdout,"-tol    epsilon          tolerance (default %E) \n",info.tol);

    exit(1);
}

static int parsar(int argc, char *argv[], pMesh mesh_distance, pMesh mesh_omega0) {
  int i, k;

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {

        case 'h':
          usage(argv[0]);
          break;

        case 'd':
          if ( !strcmp(argv[i],"-dref") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              mesh_omega0->nref1 = atoi(argv[i]);
            else
              --i;
            for (k=0; k<mesh_omega0->nref1; k++){
              ++i;
              if ( i < argc && isdigit(argv[i][0]))
                mesh_omega0->ref1[k] = atoi(argv[i]);
            }
          }
          break;

        case 'e':
          if ( !strcmp(argv[i],"-elref") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              mesh_omega0->nref2 = atoi(argv[i]);
            else
              --i;
            for (k=0; k<mesh_omega0->nref2; k++){
              ++i;
              if ( i < argc && isdigit(argv[i][0]))
                mesh_omega0->ref2[k] = atoi(argv[i]);
            }
          }
          break;

        case 'b':
          if ( !strcmp(argv[i],"-bref") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              mesh_omega0->nref3 = atoi(argv[i]);
            else
              --i;
            for (k=0; k<mesh_omega0->nref3; k++){
              ++i;
              if ( i < argc && isdigit(argv[i][0]))
                mesh_omega0->ref3[k] = atoi(argv[i]);
            }
          }
          break;

        case 'n':
          if ( !strcmp(argv[i],"-nit") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.nit = atoi(argv[i]);
            else
              --i;
          }
          break;

        case 't':
          if ( !strcmp(argv[i],"-tol") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.tol = strtod(argv[i],NULL);
            else
              --i;
          }
          break;

        default:
          fprintf(stderr,"  Unrecognized option %s\n",argv[i]);
          usage(argv[0]);
      }
    }
    else {
      if ( mesh_distance->name == NULL ) {
        mesh_distance->name = argv[i];
      }
      else if ( mesh_omega0->name == NULL ) {
        mesh_omega0->name = argv[i];
      }
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */

  if ( mesh_distance->name == NULL ) {
    fprintf(stdout,"  -- MESH DISTANCE BASENAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh_distance->name);
  }

  if ( mesh_omega0->name == NULL ) {
    fprintf(stdout,"  -- MESH OMEGA 0 BASENAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh_omega0->name);
  }

  return 1;
}

static void setfunc(int dim) {
    if ( dim == 2 ) {
        matA_P1   =  matA_P1_2d;
        hashelt   =  hashelt_2d;
        newBucket =  newBucket_2d;
        buckin    =  buckin_2d;
        locelt    =  locelt_2d;
        intpp1    =  intpp1_2d;
        evalderFunctional  = evalderFunctional2D;
        evalFunctional  = evalFunctional2D;
    }
    else {
        matA_P1   =  matA_P1_3d;
        hashelt   =  hashelt_3d;
        newBucket =  newBucket_3d;
        buckin    =  buckin_3d;
        locelt    =  locelt_3d;
        intpp1    =  intpp1_3d;
        evalderFunctional  = evalderFunctional3D;
        evalFunctional  = evalFunctional3D;
    }
}

FILE *out;

int main(int argc, char **argv) {

    Instance instance;
    int      ier;
    char     stim[32];
    pBucket  tmpBucket;
    fprintf(stdout,"  -- ELASTIC MORPHING, Release %s (%s) \n",LS_VER,LS_REL);
    fprintf(stdout,"     %s\n",LS_CPY);
    fprintf(stdout,"    %s\n",COMPIL);

    out = fopen("funct.data","a+");

    /* trap exceptions */
    signal(SIGABRT,excfun);
    signal(SIGFPE,excfun);
    signal(SIGILL,excfun);
    signal(SIGSEGV,excfun);
    signal(SIGTERM,excfun);
    signal(SIGINT,excfun);
    signal(SIGBUS,excfun);

    /* default values */

    memset(&instance.mesh_omega0,0,sizeof(Mesh));
    memset(&instance.mesh_distance,0,sizeof(Mesh));
    memset(&instance.sol_distance,0,sizeof(Sol));
    memset(&instance.sol_omega0,0,sizeof(Sol));

    info.imprim = -99;

    /* set default parameters */
    instance.mesh_omega0.nref1 = 1;
    instance.mesh_omega0.nref2 = 1;
    instance.mesh_omega0.nref3 = 1;
    instance.mesh_omega0.ref1[0] = MESH_DREF;
    instance.mesh_omega0.ref2[0] = MESH_ELREF;
    instance.mesh_omega0.ref3[0] = MESH_BREF;
    info.nit = DEFAULT_NIT;
    info.tol = DEFAULT_TOL;

    /* command line */
    if ( !parsar(argc, argv, &instance.mesh_distance, &instance.mesh_omega0) )  return 1;

    /* load data */
    if ( info.imprim )   fprintf(stdout,"\n  -- INPUT DATA\n");
    /* load meshes */
    if ( !loadMesh(&instance.mesh_distance) )  return 1;
    if ( !loadMesh(&instance.mesh_omega0) )  return 1;
     /*load signed distance function */
    instance.sol_distance.name = instance.mesh_distance.name;
    ier = loadSol(&instance.sol_distance);
    if ( !ier  )  return 1;
    if ( instance.sol_distance.np != instance.mesh_distance.np ) {
        fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. PROVIDE A VALID .sol FILE \n");
        return 0;
    }
    setfunc(instance.mesh_distance.dim);
    if (instance.mesh_distance.ne) {
        instance.mesh_distance.adja = (int*)calloc(4*instance.mesh_distance.ne+5,sizeof(int));
        assert(instance.mesh_distance.adja);
    }
    else if (instance.mesh_distance.nt) {
    	instance.mesh_distance.adja = (int*)calloc(3*instance.mesh_distance.nt+5,sizeof(int));
        assert(instance.mesh_distance.adja);
    }
    instance.sol_omega0.name = instance.mesh_omega0.name;
    instance.sol_omega0.np  = instance.mesh_omega0.np;
    instance.sol_omega0.mat = (Mat*)calloc(LS_MAT,sizeof(Mat));
    instance.sol_omega0.nit = LS_MAXIT;
    instance.sol_omega0.dim = instance.mesh_omega0.dim;
    instance.sol_omega0.ver = instance.mesh_omega0.ver;
    instance.sol_omega0.cl  = (Cl*)calloc(LS_CL,sizeof(Cl));
    instance.sol_omega0.err = LS_RES;
    instance.sol_omega0.u   = (double*)calloc(instance.sol_omega0.dim*(instance.sol_omega0.np),sizeof(double));
    assert(instance.sol_omega0.u);
    instance.sol_omega0.p   = (double*)calloc(instance.sol_omega0.dim*(instance.sol_omega0.np),sizeof(double));
    assert(instance.sol_omega0.p);
    instance.sol_omega0.d   = (double*)calloc(instance.sol_omega0.np,sizeof(double));
    assert(instance.sol_omega0.d);

    fprintf(stdout,"  -- DATA READING COMPLETED. \n");

    fprintf(stdout,"\n  %s\n   MODULE MORPHING : %s (%s)\n  %s\n",LS_STR,LS_VER,LS_REL,LS_STR);
    if ( info.imprim )   fprintf(stdout,"  -- PHASE 1 : INITIALIZATION \n");
    /* create bucket data */
    if ( !hashelt(&instance.mesh_distance) )  return 1;
    tmpBucket = newBucket(&instance.mesh_distance,BUCKSIZ);
    instance.bucket = *tmpBucket;
    if ( info.imprim ) fprintf(stdout,"  -- PHASE 1 COMPLETED. \n\n");

    /* minimization */
    if ( info.imprim ) fprintf(stdout,"  -- PHASE 2 : MINIMIZATION  \n \n");
    if ( !minimization(&instance) )  return 1;
    if ( info.imprim ) fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);

    /* free mem */
    free(instance.mesh_omega0.tetra);
    free(instance.mesh_omega0.point);
    free(instance.mesh_omega0.tria);
    free(instance.mesh_omega0.edge);
    free(instance.mesh_distance.tetra);
    free(instance.mesh_distance.point);
    free(instance.mesh_distance.tria);
    free(instance.mesh_distance.adja);
    free(instance.mesh_distance.edge);
    free(instance.sol_distance.u);
    free(instance.sol_distance.valp1);
    free(instance.sol_omega0.u);
    free(instance.sol_omega0.d);
    free(instance.sol_omega0.mat);
    free(instance.sol_omega0.cl);
    free(instance.sol_omega0.p);
    free(instance.bucket.head);
    free(instance.bucket.link);
    free(tmpBucket);

    fprintf(stdout,"\n \n   END OF MODULE MORPHING.\n \n");


    return 0;
}
