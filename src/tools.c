#include "morphing.h"

extern unsigned char inxt2[5];
extern Info  info;


int  computedistance(pInstance instance){
  int             k, icase, ret, ielem;
  double          cb[4], dist, *a, q[3];
  
  
  memset(instance->sol_omega0.d, 0, instance->sol_omega0.np*sizeof(double));
  
  for (k = 1; k <= instance->mesh_omega0.np; k++) {
    a = &instance->mesh_omega0.point[k].c[0];
    q[0] = a[0];
    q[1] = a[1];
    q[2] = a[2];
    
    icase  = buckin(&instance->mesh_distance,&instance->bucket,q);
    ielem  = locelt(&instance->mesh_distance,icase,q,cb);
    if (ielem) {
      if ( instance->mesh_distance.dim==3)
        ret = intpp1(&instance->sol_distance,instance->mesh_distance.tetra[ielem].v,&dist,ielem,cb);
      else
        ret = intpp1(&instance->sol_distance,instance->mesh_distance.tria[ielem].v,&dist,ielem,cb);
      if (!ret)  exit(1);
      instance->sol_omega0.d[k-1] = dist;
    }
  }
  
  return 1;
}


double  errL2_2d(pInstance instance){
  int             k, icase, ret, ielem, ref;
  double          cb[4], dist, *a,*b, q[2], errL2, meas, len;
  pEdge           pa;
  
  len=0;
  errL2 = 0;
  meas = 0;
  
  for (k = 1; k <= instance->mesh_omega0.na; k++) {
    pa = &instance->mesh_omega0.edge[k];
    ref = getRefm(&instance->mesh_omega0,pa->ref);
    if (ref) {
      a = &instance->mesh_omega0.point[pa->v[0]].c[0];
      b = &instance->mesh_omega0.point[pa->v[1]].c[0];
      
      len = sqrt( (b[0]- a[0])*(b[0]- a[0]) + (b[1]- a[1])*(b[1]- a[1]) );
      
      q[0] = (a[0]+b[0])/2;
      q[1] = (a[1]+b[1])/2;
      
      
      icase  = buckin(&instance->mesh_distance, &instance->bucket, q);
      ielem  = locelt(&instance->mesh_distance, icase, q, cb);
      
      if (ielem) {
        ret = intpp1(&instance->sol_distance,instance->mesh_distance.tria[ielem].v,&dist,ielem,cb);
        if (!ret)  exit(1);
      }
      
      errL2 += len*dist*dist;
      meas += len;
      
    }
  }
  
  errL2 /= meas;
  errL2 = sqrt(errL2);
  
  
  return errL2;
}

double errL2_3d(pInstance instance){
  
  pTria  pt;
  int    k, icase, ret, ielem, i, ref;
  double cb[4], *a, *b, *c, q[3], area, dist, n[3], errL2, meas;
  
  errL2 = 0;
  meas = 0;
  
  
  for (k = 1; k <= instance->mesh_omega0.nt; k++) {
    pt = &instance->mesh_omega0.tria[k];
    ref = getRefm(&instance->mesh_omega0,pt->ref);
    if (ref) {
      a = &instance->mesh_omega0.point[pt->v[0]].c[0];
      b = &instance->mesh_omega0.point[pt->v[1]].c[0];
      c = &instance->mesh_omega0.point[pt->v[2]].c[0];
      
      area = area_3d(a, b, c, n);
      
      for (i = 0; i < 3; i++) q[i] = (a[i]+b[i]+c[i])/3;
      
      icase  = buckin(&instance->mesh_distance, &instance->bucket, q);
      ielem  = locelt(&instance->mesh_distance, icase, q, cb);
      
      if (ielem) {
        ret = intpp1(&instance->sol_distance, instance->mesh_distance.tetra[ielem].v, &dist, ielem, cb);
        if (!ret)  exit(1);
      }
      errL2 += area*dist*dist;
      meas += area;
      
    }
  }
  errL2 *= 1./meas;
  
  return sqrt(errL2);
}


/* translation of centre c[3] and scaling of factor R */
int initMesh_3d(pMesh mesh, double c[3], double R){
  pPoint   ppt;
  int      k,i,m;
  double   min[3], max[3],o[3],Ray[3], rMax;
  
  /* compute center and ray mesh */
  for (i=0; i<3; i++) {
    min[i] =  FLT_MAX;
    max[i] = -FLT_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<3; i++) {
      if ( ppt->c[i] > max[i] )    max[i] = ppt->c[i];
      if ( ppt->c[i] < min[i] )    min[i] = ppt->c[i];
    }
  }
  for (m=0; m<3; m++) {
    o[m] = 0.5*( max[m] + min[m] );
    Ray[m] = max[m] - o[m];
  }
  
  rMax = Ray[0];
  rMax = LS_MAX(rMax,Ray[1]);
  rMax = LS_MAX(rMax,Ray[2]);
  
  
  /* translation + omotetie */
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (m=0; m<3; m++) {
      ppt->c[m] = (R /rMax) * (ppt->c[m]-o[m]) + c[m];
    }
  }
  
  return 1;
}


/* translation of centre c[3] and scaling of factor R */
int initMesh_2d(pMesh mesh, double c[2], double R){
  pPoint   ppt;
  int      k,i,m;
  double   min[2], max[2],o[2],Ray[2];
  
  /* compute center and ray mesh */
  for (i=0; i<2; i++) {
    min[i] =  FLT_MAX;
    max[i] = -FLT_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<2; i++) {
      if ( ppt->c[i] > max[i] )    max[i] = ppt->c[i];
      if ( ppt->c[i] < min[i] )    min[i] = ppt->c[i];
    }
  }
  for (m=0; m<2; m++) {
    o[m] = 0.5*( max[m] + min[m] );
    Ray[m] = max[m] - o[m];
  }
  
  
  /* translation + omotetie */
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (m=0; m<2; m++) {
      ppt->c[m] = (R /Ray[m]) * (ppt->c[m]-o[m]) + c[m];
    }
  }
  return 1;
}

/* return squared distance from point pa to segment (p1,p2);
 pa->tag set to 2 if distance is realized by p1 or p2 */
static double distpt_2d(pPoint p1,pPoint p2,pPoint pa,int *proj) {
  double   a,b,c,d,dd,ux,uy,vx,vy,wx,wy,xp,yp,lambda;
  
  
  *proj = 1;
  
  a = p1->c[1] - p2->c[1];
  b = p2->c[0] - p1->c[0];
  c = -b*p1->c[1] - a*p1->c[0];
  d = INIVAL_2d;
  
  dd = a*a + b*b;
  if ( dd < EPS1 ) {
    d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    *proj =2;
    return(d);
  }
  
  lambda = (pa->c[0] - p2->c[0])*(-b) + (pa->c[1] - p2->c[1])*a;
  dd = 1.0 / dd;
  lambda *=dd;
  xp = p2->c[0] + lambda*(-b);
  yp = p2->c[1] + lambda*a;
  
  ux = xp - p1->c[0];
  uy = yp - p1->c[1];
  vx = xp - p2->c[0];
  vy = yp - p2->c[1];
  wx = p2->c[0] - p1->c[0];
  wy = p2->c[1] - p1->c[1];
  
  if ( fabs(b) < EPS1 ) {
    if ( uy*wy <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
      *proj = 2;
    }
    else if ( vy*wy <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
      *proj = 2;
    }
  }
  else {
    if ( ux*wx <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
      *proj = 2;
    }
    else if ( vx*wx <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
      *proj = 2;
    }
  }
  
  return(d);
}


/* Compute Hausdorff distance between two set of edges already appearing in mesh */
double hausdorff_2d(pMesh mesh1, pMesh mesh2){
  pEdge pa,pat;
  pPoint p0,p1,pmil,p2,p3;
  Point mil;
  double rho1,rho2,haus,d, d0,d1,dmil;
  int k,j,proj,nac1,nac2;
  
  rho1 = 0.0;
  rho2 = 0.0;
  pmil = &mil;
  
  nac1 = 0;
  nac2 = 0;
  
  for(k=1; k<=mesh1->na;k++){
    pa = &mesh1->edge[k];
    pa->ref = 0;
  }
  
  for(k=1; k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    pa->ref = 0;
  }
  
  /* Set active edges for mesh 1 (discard bounding box) */
  for(k = 1; k<=mesh1->na;k++){
    pa = &mesh1->edge[k];
    p0 = &mesh1->point[pa->v[0]];
    p1 = &mesh1->point[pa->v[1]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)|| (p0->c[1]<0.01)||(p0->c[1]>0.99))
      continue;
    pa->ref  =1;
    nac1++;
  }
  
  /* Set active edges for mesh 2 */
  for(k = 1; k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    p0 = &mesh2->point[pa->v[0]];
    p1 = &mesh2->point[pa->v[1]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)|| (p0->c[1]<0.01)||(p0->c[1]>0.99))
      continue;
    
    pa->ref  =1;
    nac2++;
  }
  
  printf("Number of active edges %d %d \n", nac1, nac2);
  
  /* Compute rho(\Gamma_1,\Gamma_2)*/
  for(k=1;k<=mesh1->na;k++){
    pa = &mesh1->edge[k];
    
    if(!pa->ref) continue;
    
    p0 = &mesh1->point[pa->v[0]];
    p1 = &mesh1->point[pa->v[1]];
    pmil->c[0] = 0.5*(p0->c[0]+p1->c[0]);
    pmil->c[1] = 0.5*(p0->c[1]+p1->c[1]);
    
    d0 = 10.0;
    d1 = 10.0;
    dmil = 10.0;
    
    for(j=1;j<=mesh2->na;j++){
      pat = &mesh2->edge[j];
      p2 = &mesh2->point[pat->v[0]];
      p3 = &mesh2->point[pat->v[1]];
      
      d = distpt_2d(p2,p3,p0,&proj);
      d0 = LS_MIN(d0,d);
      
      d = distpt_2d(p2,p3,p1,&proj);
      d1 = LS_MIN(d1,d);
      
      d = distpt_2d(p2,p3,pmil,&proj);
      dmil = LS_MIN(dmil,d);
      
    }
    rho1 = LS_MAX(rho1,d0);
    rho1 = LS_MAX(rho1,d1);
    rho1 = LS_MAX(rho1,dmil);
  }
  
  /* Compute rho(\Gamma_2,\Gamma_1) */
  for(k=1;k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    
    if(!pa->ref) continue;
    
    p0 = &mesh2->point[pa->v[0]];
    p1 = &mesh2->point[pa->v[1]];
    pmil->c[0] = 0.5*(p0->c[0]+p1->c[0]);
    pmil->c[1] = 0.5*(p0->c[1]+p1->c[1]);
    
    d0 = 10.0;
    d1 = 10.0;
    dmil = 10.0;
    
    for(j=1;j<=mesh1->na;j++){
      pat = &mesh1->edge[j];
      
      if(!pat->ref) continue;
      
      p2 = &mesh1->point[pat->v[0]];
      p3 = &mesh1->point[pat->v[1]];
      
      d = distpt_2d(p2,p3,p0,&proj);
      d0 = LS_MIN(d0,d);
      
      d = distpt_2d(p2,p3,p1,&proj);
      d1 = LS_MIN(d1,d);
      
      d = distpt_2d(p2,p3,pmil,&proj);
      dmil = LS_MIN(dmil,d);
    }
    rho2 = LS_MAX(rho1,d0);
    rho2 = LS_MAX(rho1,d1);
    rho2 = LS_MAX(rho1,dmil);
  }
  
  haus = LS_MAX(rho1,rho2);
  haus = sqrt(haus);
  
  return(haus);
}


/* compute L^2 error when integrating with the real distance */
/* mesh2 : box , mesh1 : template mesh */
double errdist_2d(pMesh mesh, pMesh mesh2){
  int k,l,proj,i0,i1;
  pEdge pe,pa;
  pPoint p0,p1,p2,p3,pmil;
  Point mil;
  double errL2,dist,len=0;
  
  errL2     = 0.0;
  pmil = &mil;
  
  /* Set active edges for mesh 2 (discard bounding box) */
  for(k = 1; k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    p0 = &mesh2->point[pa->v[0]];
    p1 = &mesh2->point[pa->v[1]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)|| (p0->c[1]<0.01)||(p0->c[1]>0.99))
      continue;
    pa->ref  =1;
  }
  
  
  /* Computation of L^2 errors */
  for(k=1;k<=mesh->na;k++){
    pa = &mesh->edge[k];
    if ( pa->ref != mesh->refb) { // mesh->refb is the Dirichlet reference
      i0 = pa->v[0];
      i1 = pa->v[1];
      
      p0 = &mesh->point[i0];
      p1 = &mesh->point[i1];
      pmil->c[0] = 0.5*(p0->c[0]+p1->c[0]);
      pmil->c[1] = 0.5*(p0->c[1]+p1->c[1]);
      
      dist = FLT_MAX;
      
      for(l=1;l<=mesh2->na;l++){
        pe  = &mesh2->edge[l];
        if((pe->ref)==1) {
          p2  = &mesh2->point[pe->v[0]];
          p3  = &mesh2->point[pe->v[1]];
          
          dist = LS_MIN(dist,distpt_2d(p2,p3,pmil,&proj));
        }
      }
      
      errL2 += len*dist;
    }
  }
  errL2 = sqrt(errL2);
  printf("L2 error  : \t %e\n",errL2);
  return(errL2);
}


/* Compute Hausdorff distance between two set of triangles discarding bounding box and interior Dirichlet boundary */
/* mesh1 : target mesh, mesh2 : template mesh */
double hausdorff_3d(pMesh mesh1, pMesh mesh2){
  pTria   pt,pt1;
  pPoint  p0,p1,pmil,p2,p3,p4,p5;
  Point   mil;
  double  rho2,haus,d, d0,d1,d2,dmil;
  double rho1;
  int     k,i,j,nac1,nac2;
  char    proj;
  
  rho1 = 0.0;
  rho2 = 0.0;
  pmil = &mil;
  
  nac1 = 0;
  nac2 = 0;
  
  for(k=1; k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    if ( pt->ref != mesh2->refb)
      pt->ref = 0;
  }
  
  for(k=1; k<=mesh2->nt;k++){
    pt1 = &mesh2->tria[k];
    if ( pt1->ref != mesh2->refb)
      pt1->ref = 0;
  }
  
  /* Set active triangles for mesh 1 (discard bounding box) */
  for(k = 1; k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    p0 = &mesh1->point[pt->v[0]];
    p1 = &mesh1->point[pt->v[1]];
    p2 = &mesh1->point[pt->v[2]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)||(p0->c[1]<0.01)||(p0->c[1]>0.99)||(p0->c[2]<0.01)||(p0->c[2]>0.99)) continue;
    pt->ref  =5;
    nac1++;
  }
  
  /* Set active triangles for mesh 2 (discarding Dirichlet) */
  for(k = 1; k<=mesh2->nt;k++){
    pt = &mesh2->tria[k];
    if ( pt->ref != mesh2->refb) {
      pt->ref  =5;
      nac2++;
    }
  }
  
  printf("Number of active triangles %d %d \n", nac1, nac2);
  
  /* Compute rho(\Gamma_1,\Gamma_2)*/
  for(k=1;k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    if(pt->ref ==5){
      
      p0 = &mesh1->point[pt->v[0]];
      p1 = &mesh1->point[pt->v[1]];
      p2 = &mesh1->point[pt->v[2]];
      
      for (i=0; i<3;i++) pmil->c[i] = (p0->c[i]+p1->c[i]+p2->c[i])/3;
      
      d0 = 10.0;
      d1 = 10.0;
      d2 = 10.0;
      dmil = 10.0;
      
      
      for(j=1;j<=mesh2->nt;j++){
        pt1 = &mesh2->tria[j];
        if ( pt1->ref ==5 ) {
          p3 = &mesh2->point[pt1->v[0]];
          p4 = &mesh2->point[pt1->v[1]];
          p5 = &mesh2->point[pt1->v[2]];
          
          d = distpt_3d(p3,p4,p5,p0,&proj);
          d0 = LS_MIN(d0,d);
          
          d = distpt_3d(p3,p4,p5,p1,&proj);
          d1 = LS_MIN(d1,d);
          
          d = distpt_3d(p3,p4,p5,p2,&proj);
          d2 = LS_MIN(d2,d);
          
          d = distpt_3d(p3,p4,p5,pmil,&proj);
          dmil = LS_MIN(dmil,d);
          
        }
      }
    }
    rho1 = LS_MAX(rho1,d0);
    rho1 = LS_MAX(rho1,d1);
    rho1 = LS_MAX(rho1,d2);
  }
  printf("rho1 %e  \n", rho1);
  
  /* Compute rho(\Gamma_1,\Gamma_2)*/
  for(k=1;k<=mesh2->nt;k++){
    pt = &mesh2->tria[k];
    if(pt->ref ==5){
      
      p0 = &mesh2->point[pt->v[0]];
      p1 = &mesh2->point[pt->v[1]];
      p2 = &mesh2->point[pt->v[2]];
      
      for (i=0; i<3;i++) pmil->c[i] = (p0->c[i]+p1->c[i]+p2->c[i])/3;
      
      d0 = 10.0;
      d1 = 10.0;
      d2 = 10.0;
      dmil = 10.0;
      
      
      for(j=1;j<=mesh1->nt;j++){
        pt1 = &mesh1->tria[j];
        p3 = &mesh1->point[pt1->v[0]];
        p4 = &mesh1->point[pt1->v[1]];
        p5 = &mesh1->point[pt1->v[2]];
        
        d = distpt_3d(p3,p4,p5,p0,&proj);
        d0 = LS_MIN(d0,d);
        
        d = distpt_3d(p3,p4,p5,p1,&proj);
        d1 = LS_MIN(d1,d);
        
        d = distpt_3d(p3,p4,p5,p2,&proj);
        d2 = LS_MIN(d2,d);
        
        d = distpt_3d(p3,p4,p5,pmil,&proj);
        dmil = LS_MIN(dmil,d);
        
        
      }
    }
    rho2 = LS_MAX(rho1,d0);
    rho2 = LS_MAX(rho1,d1);
    rho2 = LS_MAX(rho1,d2);
    rho2 = LS_MAX(rho1,dmil);
  }
  
  haus = LS_MAX(rho1,rho2);
  haus = sqrt(haus);
  
  return(haus);
}


/* Compute squared distance from pq to tria p0p1p2,
 proj = 1: projection onto face, =2: distance to vertex or edge */
double distpt_3d(pPoint p0,pPoint p1,pPoint p2,pPoint pq,char *proj) {
  Point pointPlan0, pointPlan1, pointPlan2, pointPlana;
  double lx1, ly1, lz1, lx2, ly2, lz2, lxq, lyq, lzq, longp0p1, cosAlpha, sinAlpha, aa, bb, ab, ll, l;
  double p0XTemp, p0YTemp, p0ZTemp, p1XTemp, p1YTemp, p1ZTemp, p2XTemp, p2YTemp, p2ZTemp, pqXTemp, pqYTemp, pqZTemp; // pour stocker temporairement
  double m11, m12, m13, m21, m22, m23, m31, m32, m33, d01, d12, d02, dTmp;
  double zone01, zone02, zone12;
  
  *proj = 1;
  lx1 = p1->c[0] - p0->c[0];
  ly1 = p1->c[1] - p0->c[1];
  lz1 = p1->c[2] - p0->c[2];
  
  lx2 = p2->c[0] - p0->c[0];
  ly2 = p2->c[1] - p0->c[1];
  lz2 = p2->c[2] - p0->c[2];
  
  lxq = pq->c[0] - p0->c[0];
  lyq = pq->c[1] - p0->c[1];
  lzq = pq->c[2] - p0->c[2];
  
  longp0p1 = sqrt(lx1*lx1 + ly1*ly1 + lz1*lz1);
  
  cosAlpha = lz1/longp0p1;
  sinAlpha = sqrt(1.0-cosAlpha*cosAlpha);
  
  p0XTemp = 0.0;
  p0YTemp = 0.0;
  p0ZTemp = 0.0;
  p1XTemp = lx1;
  p1YTemp = ly1;
  p1ZTemp = lz1;
  p2XTemp = lx2;
  p2YTemp = ly2;
  p2ZTemp = lz2;
  pqXTemp = lxq;
  pqYTemp = lyq;
  pqZTemp = lzq;
  
  /* Apply rotation with angle alpha */
  
  aa = lx1*lx1;
  bb = ly1*ly1;
  ab = lx1*ly1;
  ll = aa + bb;
  l = sqrt(ll);
  
  if(ll<EPS1){
    if(p1ZTemp <= 0.0){
      p1ZTemp = -p1ZTemp;
      p2ZTemp = -p2ZTemp;
      pqZTemp = -pqZTemp;
    }
  }
  
  else{
    m11 = (aa*cosAlpha + bb)/ll;
    m12 = (ab*cosAlpha - ab )/ll;
    m13 = -lx1*sinAlpha/l;
    
    m21 = (ab * cosAlpha - ab)/ll;
    m22 = (bb*cosAlpha + aa)/ll;
    m23 = - ly1 * sinAlpha/l;
    
    m31 = (lx1* sinAlpha)/l;
    m32 = (ly1 * sinAlpha)/l;
    m33 = cosAlpha;
    
    p1XTemp = m11*lx1 + m12*ly1 + m13 * lz1;
    p1YTemp = m21*lx1 + m22*ly1 + m23 * lz1;
    p1ZTemp = longp0p1;
    
    p2XTemp = m11*lx2 + m12*ly2 + m13 * lz2;
    p2YTemp = m21*lx2 + m22*ly2 + m23 * lz2;
    p2ZTemp = m31*lx2 + m32*ly2 + m33 * lz2;
    
    pqXTemp = m11*lxq + m12*lyq + m13 * lzq;
    pqYTemp = m21*lxq + m22*lyq + m23 * lzq;
    pqZTemp = m31*lxq + m32*lyq + m33 * lzq;
  }
  
  /* Apply second rotation to put point p2 in plane (yz), in half plane y > 0*/
  assert((p2XTemp* p2XTemp + p2YTemp* p2YTemp) > 0.0);  // au cas ou...
  cosAlpha = p2YTemp/(sqrt(p2XTemp* p2XTemp + p2YTemp* p2YTemp));
  sinAlpha = sqrt(1.0-cosAlpha * cosAlpha);
  if(p2XTemp <=0.0) {sinAlpha = -sinAlpha; }
  
  lx2 = 0.0;
  ly2 = sinAlpha * p2XTemp + cosAlpha * p2YTemp;
  lz2 = p2ZTemp;
  
  lxq = cosAlpha * pqXTemp - sinAlpha * pqYTemp;
  lyq = sinAlpha * pqXTemp + cosAlpha * pqYTemp;
  lzq = pqZTemp;
  
  zone01 = lyq;
  zone02 = ly2*lzq - lz2*lyq;
  zone12 = lyq*(lz2-longp0p1) - ly2*(lzq - longp0p1);
  
  if((zone01>=0.0)&&(zone02 >=0.0)&&(zone12>=0.0))
  {
    return (lxq*lxq);
  }
  
  else
  {
    pointPlan0.c[0] = 0.0;
    pointPlan0.c[1] = 0.0;
    
    pointPlan1.c[0] = 0.0;
    pointPlan1.c[1] = longp0p1;
    
    pointPlan2.c[0] = ly2;
    pointPlan2.c[1] = lz2;
    
    pointPlana.c[0] = lyq;
    pointPlana.c[1] = lzq;
    
    d01 = distpt_23d(&pointPlan0, &pointPlan1, &pointPlana);
    d02 = distpt_23d(&pointPlan0, &pointPlan2, &pointPlana);
    d12 = distpt_23d(&pointPlan1, &pointPlan2, &pointPlana);
    
    dTmp = LS_MIN(d01, LS_MIN(d02, d12));
    *proj = 2;
    
    return(dTmp + lxq*lxq);
  }
  
}

/* return distance from point pa to segment (p1p2)
 vertex tag set to 2 if distance is realized by p1 or p2 */
double distpt_23d(pPoint p1,pPoint p2,pPoint pa) {
  double   a,b,c,d,dd,ux,uy,vx,vy,wx,wy,xp,yp;
  
  a = p1->c[1] - p2->c[1];
  b = p2->c[0] - p1->c[0];
  c = -b*p1->c[1] - a*p1->c[0];
  d = INIVAL_3d;
  
  dd = a*a + b*b;
  if ( dd < EPS1 ) {
    d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    return(d);
  }
  xp =  b*b * pa->c[0] - a*b * pa->c[1] - a*c;
  yp = -a*b * pa->c[0] + a*a * pa->c[1] - b*c;
  dd = 1.0 / dd;
  xp *= dd;
  yp *= dd;
  
  ux = xp - p1->c[0];
  uy = yp - p1->c[1];
  vx = xp - p2->c[0];
  vy = yp - p2->c[1];
  wx = p2->c[0] - p1->c[0];
  wy = p2->c[1] - p1->c[1];
  
  if ( fabs(b) < EPS1 ) {
    if ( uy*wy <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    }
    else if ( vy*wy <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
    }
  }
  else {
    if ( ux*wx <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    }
    else if ( vx*wx <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
    }
  }
  
  return(d);
}


/* compute L^2 error when integrating with the real distance */
/* mesh2 : box , mesh1 : template mesh */
double errdist_3d(pMesh mesh1, pMesh mesh2){
  int k,l,nac1,nac2;
  pTria pt,pt1;
  pPoint p0,p1,p2,p3,p4,p5,pmil,pmil1,pmil2,pmil3;
  Point mil;
  double errL2,dist,area,n[3],d1,d2,d3,meas;
  char proj;
  
  nac1 = 0.0;
  nac2 = 0.0;
  meas = 0.0;
  
  errL2     = 0.0;
  pmil = &mil;
  pmil1 = &mil;
  pmil2 = &mil;
  pmil3 = &mil;
  
  
  /* Set active triangles for mesh 1 (discard bounding box) */
  for(k = 1; k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    p0 = &mesh1->point[pt->v[0]];
    p1 = &mesh1->point[pt->v[1]];
    p2 = &mesh1->point[pt->v[2]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)||(p0->c[1]<0.01)||(p0->c[1]>0.99)||(p0->c[2]<0.01)||(p0->c[2]>0.99)) continue;
    pt->ref  =5;
    nac1++;
  }
  
  /* Set active triangles for mesh 2 (discarding Dirichlet) */
  for(k = 1; k<=mesh2->nt;k++){
    pt = &mesh2->tria[k];
    if ( pt->ref != mesh2->refb) {
      pt->ref =5;
      nac2++;
    }
  }
  
  printf("Number of active triangles %d %d \n", nac1, nac2);
  
  /* Computation of L^2 errors */
  for(k=1;k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    
    if ( pt->ref ==5) { // mesh->ref is the Dirichlet reference
      
      printf("Tria k  %d\n", k);
      
      p0 = &mesh1->point[pt->v[0]];
      p1 = &mesh1->point[pt->v[1]];
      p2 = &mesh1->point[pt->v[2]];
      
      
      pmil->c[0] = (0.5*p0->c[0]+0.5*p1->c[0]);
      pmil->c[1] = (0.5*p0->c[1]+0.5*p1->c[1]);
      pmil->c[2] = (0.5*p0->c[2]+0.5*p1->c[2]);
      
      
      pmil1->c[0] = (0.5*p0->c[0]+0.5*p2->c[0]);
      pmil1->c[1] = (0.5*p0->c[1]+0.5*p2->c[1]);
      pmil1->c[2] = (0.5*p0->c[2]+0.5*p2->c[2]);
      
      pmil2->c[0] = 0.5*(p1->c[0]+p2->c[0]);
      pmil2->c[1] = 0.5*(p1->c[1]+p2->c[1]);
      pmil2->c[2] = 0.5*(p1->c[2]+p2->c[2]);
      
      
      dist = FLT_MAX;
      d1 = FLT_MAX;
      d2 = FLT_MAX;
      d3 = FLT_MAX;
      
      for(l=1;l<=mesh2->nt;l++){
        pt1  = &mesh2->tria[l];
        if(pt1->ref ==5) {
          p3  = &mesh2->point[pt1->v[0]];
          p4  = &mesh2->point[pt1->v[1]];
          p5  = &mesh2->point[pt1->v[2]];
          
          d1 = LS_MIN(d1,distpt_3d(p3,p4,p5,pmil,&proj));
          d2 = LS_MIN(d2,distpt_3d(p3,p4,p5,pmil1,&proj));
          d3 = LS_MIN(d3,distpt_3d(p3,p4,p5,pmil2,&proj));
          //dist = LS_MIN(dist,distpt_3d(p3,p4,p5,pmil,&proj));
        }
      }
      
      area = area_3d(p3->c,p4->c,p5->c,n);
      
      errL2 += area*( d1*d1 + d2*d2 + d3*d3)/6;
      meas += area;
    }
  }
  meas = 1./meas;
  errL2 = sqrt(errL2*meas);
  printf("L2 error  : \t %e\n",errL2);
  return(errL2);
}