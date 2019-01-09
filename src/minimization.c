#include "morphing.h"

#define BUCKSIZ 64

extern FILE *out;
extern Info info;

/* retrieve references in list */
int getRefm(pMesh mesh, int ref) {
  int     i;

  for (i=0; i<mesh->nref3; i++) {
    if (mesh->ref3[i] == ref )  return(1);
  }

  return 0;
}


/* find best step */
int findstep(pInstance instance, double *step){
  double  J,J_temp;

  J = evalFunctional(instance,0);
  J_temp = evalFunctional3D(instance,step[0]);

  while (J_temp>=J) {
    step[0] *= 0.8;
    J_temp = evalFunctional(instance,step[0]);
    if (step[0] < 1e-8) {
      return 0;
    }
  }
  return 1;
}

int minimization(pInstance instance){

  int    k=0;
  double *F, step, tol;
  step = 1e8;

  F = (double*)calloc(instance->mesh_omega0.dim*instance->mesh_omega0.np, sizeof(double));
  memset(instance->sol_omega0.u, 0, instance->mesh_omega0.dim*instance->sol_omega0.np*sizeof(double));
  assert(F);

  /* compute initial distance field */
  fprintf(stdout," \n Writing initial distance field ..\n");
  if ( !computedistance(instance) )
    return 0;

  if ( !saveDistance(&instance->sol_omega0, &instance->mesh_omega0, &instance->mesh_distance,0) )
      return 0;

  /*
  if(instance->mesh_omega0.dim==3){
    fprintf(stdout, "if ok\n");
    tol=errL2_3d(instance);
  }
  else{
    fprintf(stdout, "if pas ok\n");
    tol=errL2_3d(instance);
  }
  */

  /* gradient descent iterations */
  while(k<=info.nit ){// && tol > info.tol){
    fprintf(stdout, "Iteration: %d / %d\n", k, info.nit);
    /* compute descent direction */
    if ( !evalderFunctional(instance,F) )  return 0;
    /* find best step */
    if ( !findstep(instance,&step) ) step = 1e3;
    /* advect mesh */
    if ( !moveMesh(&instance->mesh_omega0,&instance->sol_omega0,step) )  return 0;

    /*
    if(k%20==0) { step = 1e7;
      if(instance->mesh_omega0.dim==3) tol=errL2_3d(instance);
      else tol=errL2_3d(instance);
    }
    */
    k++;

  }

  if(k==info.nit+1) fprintf(stdout," \n Max number of iterations reached it = %d .. \n",info.nit);
  else if (tol > info.tol) fprintf(stdout," \n Desired tolerance reached tol = %e .. \n",info.tol);

  fprintf(stdout," \n Writing data files .. \n");
  if ( !saveMesh(&instance->mesh_omega0, &instance->mesh_distance,1) )  return 0;
  if ( !computedistance(instance) )  return 0;
  if ( !saveDistance(&instance->sol_omega0, &instance->mesh_omega0, &instance->mesh_distance,1))  return 0;
  if ( !saveSol(&instance->sol_omega0, &instance->mesh_omega0, &instance->mesh_distance, 1) )  return 0;
  fprintf(stdout,"\n");

  /*free memory */
  free(F);

  return 1;
}
