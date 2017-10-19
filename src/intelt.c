#include "morphing.h"



/* P1 interpolation */
int intpp1_3d(pSol sol1,int *v,double *sp,int ip,double *cb) {
    double  *sa,*sb,*sc,*sd;
    int      i,iada,iadb,iadc,iadd;
    
    
    
    iada = (v[0]-1)*sol1->size[0] ;
    iadb = (v[1]-1)*sol1->size[0] ;
    iadc = (v[2]-1)*sol1->size[0] ;
    iadd = (v[3]-1)*sol1->size[0] ;
  
//      iada = (v[0]-1)*sol1->size[0] +1;
//      iadb = (v[1]-1)*sol1->size[0] +1;
//      iadc = (v[2]-1)*sol1->size[0] +1;
//      iadd = (v[3]-1)*sol1->size[0] +1;

    sa = &sol1->valp1[iada];
    sb = &sol1->valp1[iadb];
    sc = &sol1->valp1[iadc];
    sd = &sol1->valp1[iadd];
    
    
    for (i=0; i<sol1->size[0]; i++)
        
        sp[i] = cb[0]*sa[i] + cb[1]*sb[i] + cb[2]*sc[i] + cb[3]*sd[i];
    
    return(1);
}


/* P1 interpolation */
int intpp1_2d(pSol sol1,int *v,double *sp,int ip,double *cb) {
    double  *sa,*sb,*sc;
    int      i,iada,iadb,iadc;
  

 
//    iada = (v[0]-1)*sol1->size[0] + 1;
//    iadb = (v[1]-1)*sol1->size[0] + 1;
//    iadc = (v[2]-1)*sol1->size[0] + 1;
  
  iada = (v[0]-1)*sol1->size[0] ;
  iadb = (v[1]-1)*sol1->size[0] ;
  iadc = (v[2]-1)*sol1->size[0] ;
 
  
  
    sa = &sol1->valp1[iada];
    sb = &sol1->valp1[iadb];
    sc = &sol1->valp1[iadc];

    
    for (i=0; i<sol1->size[0]; i++)
        sp[i] = cb[0]*sa[i] + cb[1]*sb[i] + cb[2]*sc[i];

    
    return(1); 
}




