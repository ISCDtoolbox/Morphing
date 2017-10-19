#include "morphing.h"

extern Info    info;



int scaleMesh(pMesh mesh) {
  pPoint    ppt;
  double    dd,delta,min[3],max[3];
  int       i,k;
  
  /* compute bounding box */
  for (i=0; i<mesh->dim; i++) {
    min[i] =  FLT_MAX;
    max[i] = -FLT_MAX;
  }
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > max[i] )  max[i] = ppt->c[i];
      if ( ppt->c[i] < min[i] )  min[i] = ppt->c[i];
    }
  }
  delta = 0.0;
  for (i=0; i<mesh->dim; i++) {
    dd = fabs(max[i]-min[i]);
    if ( dd > delta )  delta = dd;
  }
  if ( delta < EPSD ) {
    fprintf(stdout,"  ## Unable to scale mesh\n");
    return(0);
  }
  
  /* normalize coordinates */
  /*dd = (double)PRECI / info.delta;*/
  dd = 0.9;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++)
      ppt->c[i] = (dd/delta) * (ppt->c[i] - min[i]) + 0.05;
  }
  
  /* compute new bounding box */
  for (i=0; i<mesh->dim; i++) {
    min[i] =  FLT_MAX;
    max[i] = -FLT_MAX;
  }
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > max[i] )  max[i] = ppt->c[i];
      if ( ppt->c[i] < min[i] )  min[i] = ppt->c[i];
    }
  }
fprintf(stdout,"  -- Max = %lf %lf %lf \n",max[0],max[1],max[2]);
fprintf(stdout,"  -- Min = %lf %lf %lf \n",min[0],min[1],min[2]);
fprintf(stdout,"  -- Centre = %lf %lf %lf \n",0.5*(max[0]+min[0]),0.5*(max[1]+min[1]),0.5*(max[2]+min[2]));
  
  return(1);
}

int moveMesh(pMesh mesh,pSol sol, double step) {
    pPoint    ppt;
    int       i,k;
    
    for (k=1; k<=mesh->np; k++){
        ppt = &mesh->point[k];
        for (i=0; i< mesh->dim; i++){
            ppt->c[i] = ppt->c[i] + step*sol->u[(mesh->dim)*(k-1)+i];
            sol->p[(mesh->dim)*(k-1)+i] += step*sol->u[(mesh->dim)*(k-1)+i];
        }
    }
    
    return(1);
}

