#include "morphing.h"

#define BUCKSIZ 64

extern FILE *out;
extern Info info;


double evalFunctional3D(pInstance instance, double step) {
  
    pTetra          ptt;
    int             k, icase, ret, ielem, i, ref;
    double          cb[4], vol, *a, *b, *c, *d, ag[3], bg[3], cg[3], dg[3], dist, q[3];
    double          J;

    /* evaluate J(g) */
    J = 0;
    for (k = 1; k <= instance->mesh_omega0.ne; k++) {
        ptt = &instance->mesh_omega0.tetra[k];
        ref = getRefel(&instance->mesh_omega0,ptt->ref);
        if(!ref){

        a = &instance->mesh_omega0.point[ptt->v[0]].c[0];
        b = &instance->mesh_omega0.point[ptt->v[1]].c[0];
        c = &instance->mesh_omega0.point[ptt->v[2]].c[0];
        d = &instance->mesh_omega0.point[ptt->v[3]].c[0];
        
        for (i = 0; i < 3; i++) {
            ag[i] = a[i] + step*instance->sol_omega0.u[3*(ptt->v[0]-1)+i];
            bg[i] = b[i] + step*instance->sol_omega0.u[3*(ptt->v[1]-1)+i];
            cg[i] = c[i] + step*instance->sol_omega0.u[3*(ptt->v[2]-1)+i];
            dg[i] = d[i] + step*instance->sol_omega0.u[3*(ptt->v[3]-1)+i];
        }
        vol = volume(ag, bg, cg, dg);
        if (vol <= 0) {
        	return 1;
        }
        
        /* Gauss rule with 1 node */
        for (i = 0; i < 3; i++) q[i] = (ag[i]+bg[i]+cg[i]+dg[i])/4;
        
        if (q[0] < 0 || q[0] > 1 || q[1] < 0 || q[1] > 1 || q[2] < 0 || q[2] > 1) {
        	J += vol;
        }
        else {

		    icase  = buckin(&instance->mesh_distance, &instance->bucket, q);
		    ielem  = locelt(&instance->mesh_distance, icase, q, cb);
	
	        if (ielem) {
	            ret = intpp1(&instance->sol_distance, instance->mesh_distance.tetra[ielem].v, &dist, ielem, cb);
	            if (!ret)  exit(1);
            
           
	        }
	
	        J += vol*dist;
        }
      }
    }
  
    return J;
}


int evalderFunctional3D(pInstance instance, double *F) {
    
    pTria  pt;
    int    k, icase, ret, ielem, i, ref, ref1;
    double cb[4], *a, *b, *c, q[3], area, dist, n[3], coef;

    
    coef=-1;
    memset(instance->sol_omega0.u, 0, 3*instance->sol_omega0.np*sizeof(double));
    memset(F, 0, 3*instance->sol_omega0.np*sizeof(double));
    /* Assembly elastic problem  */
    for (k = 1; k <= instance->mesh_omega0.nt; k++) {
        pt = &instance->mesh_omega0.tria[k];
      ref = getRefm(&instance->mesh_omega0,pt->ref);
      ref1 = getRefb(&instance->mesh_omega0,pt->ref);
      if ( (!ref) && (!ref1)) {
	        a = &instance->mesh_omega0.point[pt->v[0]].c[0];
	        b = &instance->mesh_omega0.point[pt->v[1]].c[0];
	        c = &instance->mesh_omega0.point[pt->v[2]].c[0];
	        
          area = area_3d(a, b, c, n);
	        
	        for (i = 0; i < 3; i++) q[i] = (a[i]+b[i]+c[i])/3;
	        
	        if (q[0] < 0 || q[0] > 1 || q[1] < 0 || q[1] > 1 || q[2] < 0 || q[2] > 1) {
        		for (i = 0; i < 3; i++) {
		            F[3*(pt->v[i]-1)+0] += coef*area*n[0]/3;
		            F[3*(pt->v[i]-1)+1] += coef*area*n[1]/3;
		            F[3*(pt->v[i]-1)+2] += coef*area*n[2]/3;
	      		}
        	}
      	 	else {
		        icase  = buckin(&instance->mesh_distance, &instance->bucket, q);
		        ielem  = locelt(&instance->mesh_distance, icase, q, cb);
		        
		        if (ielem) {
		            ret = intpp1(&instance->sol_distance, instance->mesh_distance.tetra[ielem].v, &dist, ielem, cb);
		            if (!ret)  exit(1);
		        }
		        for (i = 0; i < 3; i++) {
		            F[3*(pt->v[i]-1)+0] += coef*area*dist*n[0]/3;
		            F[3*(pt->v[i]-1)+1] += coef*area*dist*n[1]/3;
		            F[3*(pt->v[i]-1)+2] += coef*area*dist*n[2]/3;
		        }
          }
       }
    }
    
    /* compute the displacement u */
    if ( !elas(instance, F) )  return 0;
    
    return 1;
}


