#include "morphing.h"

#define BUCKSIZ 64

extern FILE *out;
extern Info info;

/* return length of edge + outer normal to segment */
static double length(double *a, double *b, double n[2]) {
  double ax, ay, dd;
  
  ax = b[0] - a[0];
  ay = b[1] - a[1];
  n[0] = ay;
  n[1] = -ax;
  dd   = sqrt(ax*ax+ay*ay);
  if (dd > EPSD) {
    n[0] *= 1./dd;
    n[1] *= 1./dd;
  }
  return dd;
}

/* triangle area */
static inline double area_2d(double *a, double *b, double *c) {
  double ux, uy, vx, vy, dd;
  
  ux = b[0] - a[0];
  uy = b[1] - a[1];
  vx = c[0] - a[0];
  vy = c[1] - a[1];
  dd = 0.5*(ux*vy-uy*vx);
  
  return dd;
}

double evalFunctional2D(pInstance instance, double step){
  
  double J;
  pTria  pt;
  int    k, icase, ret, ielem, i,ref;
  double cb[3], area, *a, *b, *c, ag[3], bg[3], cg[3], q[2], dist;
  
  /* evaluate J*/
  J = 0;
  
  for (k = 1; k <= instance->mesh_omega0.nt; k++) {
    pt = &instance->mesh_omega0.tria[k];
    ref = getRefel(&instance->mesh_omega0,pt->ref);
    if(!ref){
      a = &instance->mesh_omega0.point[pt->v[0]].c[0];
      b = &instance->mesh_omega0.point[pt->v[1]].c[0];
      c = &instance->mesh_omega0.point[pt->v[2]].c[0];
      
      for (i = 0; i < 2; i++) {
        ag[i] = a[i] + step*instance->sol_omega0.u[2*(pt->v[0]-1)+i];
        bg[i] = b[i] + step*instance->sol_omega0.u[2*(pt->v[1]-1)+i];
        cg[i] = c[i] + step*instance->sol_omega0.u[2*(pt->v[2]-1)+i];
      }
      
      area = area_2d(ag, bg, cg);
      if (area <= 0) {
        return 1;
      }
      
      /* Gauss rule with 1 node */
      q[0] = (ag[0]+bg[0]+cg[0])/3;
      q[1] = (ag[1]+bg[1]+cg[1])/3;
      
      if (q[0] < 0 || q[0] > 1 || q[1] < 0 || q[1] > 1) {
        J += area;
      }
      else {
        icase  = buckin(&instance->mesh_distance,&instance->bucket,q);
        ielem  = locelt(&instance->mesh_distance,icase,q,cb);
        if (ielem) {
          ret = intpp1(&instance->sol_distance,instance->mesh_distance.tria[ielem].v,&dist,ielem,cb);
          if (!ret)  exit(1);
          instance->sol_omega0.d[pt->v[0]] = dist;
        }
        
        J += area*dist;
      }
    }
  }
 fprintf(out,"%e \n",J);
  
  return J;
}

int evalderFunctional2D(pInstance instance, double *F){
  pEdge  pa;
  int    k, icase, ret, ielem, i, coef, ref, ref1;
  double cb[3],dist, n[2], len, q[2];
  double *a, *b;
  
  coef = 1;
  
  fprintf(stdout,"  --  Compute the descent direction \n");
  
  memset(instance->sol_omega0.u, 0, 2*instance->sol_omega0.np*sizeof(double));
  memset(F, 0, 2*instance->sol_omega0.np*sizeof(double));
  
  
  /* Assembly elastic second membre */
  for (k = 1; k <= instance->mesh_omega0.na; k++) {
    pa = &instance->mesh_omega0.edge[k];
    ref = getRefm(&instance->mesh_omega0,pa->ref);
    ref1 = getRefb(&instance->mesh_omega0,pa->ref);
    if ( (!ref) && (!ref1)) {
      a = &instance->mesh_omega0.point[pa->v[0]].c[0];
      b = &instance->mesh_omega0.point[pa->v[1]].c[0];
      
      len = length(a, b, n);
      
      q[0] = (a[0]+b[0])/2;
      q[1] = (a[1]+b[1])/2;
      
      if (q[0] < 0 || q[0] > 1 || q[1] < 0 || q[1] > 1) {
        for (i = 0; i < 2; i++) {
          F[2*(pa->v[i]-1)+0] -= coef*n[0]/2;
          F[2*(pa->v[i]-1)+1] -= coef*n[1]/2;
        }
      }
      else {
        icase  = buckin(&instance->mesh_distance, &instance->bucket, q);
        ielem  = locelt(&instance->mesh_distance, icase, q, cb);
        
        if (ielem) {
          ret = intpp1(&instance->sol_distance,instance->mesh_distance.tria[ielem].v,&dist,ielem,cb);
          if (!ret)  exit(1);
        }
        
        for (i = 0; i < 2; i++) {
          F[2*(pa->v[i]-1)+0] -= coef*len*dist*n[0]/2;
          F[2*(pa->v[i]-1)+1] -= coef*len*dist*n[1]/2;
        }
      }
    }
  }
  
  /* compute the displacement u */
  if ( !elas(instance, F) )  return 0;
  
  return 1;
}

