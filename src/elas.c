#include "morphing.h"

extern Info info;


/* find references in list */
int getRefb(pMesh mesh, int ref) {
  int     i;
  
  for (i=0; i<mesh->nref1; i++) {
    if (mesh->ref1[i] == ref )  return 1;
  }
  return 0;
}

/* find reference in list */
int getRefel(pMesh mesh, int ref) {
  int     i;
  
  for (i=0; i<mesh->nref2; i++) {
    if (mesh->ref2[i] == ref )  return 1;
  }
  return 0;
}

/* compute triangle area and unit normal in 3d */
double area_3d(double *a, double *b, double *c, double *n) {
    double    ux, uy, uz, vx, vy, vz, dd, dd1;
    
    ux = b[0] - a[0];
    uy = b[1] - a[1];
    uz = b[2] - a[2];
    
    vx = c[0] - a[0];
    vy = c[1] - a[1];
    vz = c[2] - a[2];
    
    n[0] = uy*vz - uz*vy;
    n[1] = uz*vx - ux*vz;
    n[2] = ux*vy - uy*vx;
    dd   = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if ( dd > EPSD ) {
        dd1 = 1.0 / dd;
        n[0] *= dd1;
        n[1] *= dd1;
        n[2] *= dd1;
    }
    
    return 0.5*dd;
}

/* compute volume of tetra */
double volume(double *a, double *b, double *c, double *d) {
    double  ax, ay, az, bx, by, bz, vol;
    
    ax = b[0] - a[0];
    ay = b[1] - a[1];
    az = b[2] - a[2];
    
    bx = c[0] - a[0];
    by = c[1] - a[1];
    bz = c[2] - a[2];
    
    vol = (d[0]-a[0]) * (ay*bz - az*by) + (d[1]-a[1]) * (az*bx - ax*bz) \
    + (d[2]-a[2]) * (ax*by - ay*bx);
    
    return (vol/6);
}


/* invert 3x3 non-symmetric matrix */
int invmatg(double m[9], double mi[9]) {
    double  aa, bb, cc, det, vmin, vmax, maxx;
    int     k;
    
    /* check ill-conditionned matrix */
    vmin = vmax = fabs(m[0]);
    for (k=1; k<9; k++) {
        maxx = fabs(m[k]);
        if ( maxx < vmin )  vmin = maxx;
        else if ( maxx > vmax )  vmax = maxx;
    }
    if ( vmax == 0. )  return 0;
    
    /* compute sub-dets */
    aa = m[4]*m[8] - m[5]*m[7];
    bb = m[5]*m[6] - m[3]*m[8];
    cc = m[3]*m[7] - m[4]*m[6];
    det = m[0]*aa + m[1]*bb + m[2]*cc;
    if ( fabs(det) < EPSD )  return 0;
    det = 1. / det;
    
    mi[0] = aa*det;
    mi[3] = bb*det;
    mi[6] = cc*det;
    mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
    mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
    mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
    mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
    mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
    mi[8] = (m[0]*m[4] - m[1]*m[3])*det;
    
    return 1;
}


static int setTGV_3d(pMesh mesh, pCsr A) {
	
	pTria    ptt;
	pTetra   pt;
    int      ig, i, k, ref;
    
    for (k = 1; k <= mesh->ne; k++) {
	    pt = &mesh->tetra[k];
      ref =getRefel(mesh,pt->ref);
	    if (ref) {
		    for (i = 0; i < 4; i++) {
		    	ig = pt->v[i];
		        csrSet(A, 3*(ig-1)+0, 3*(ig-1)+0, LS_TGV);
		        csrSet(A, 3*(ig-1)+1, 3*(ig-1)+1, LS_TGV);
		        csrSet(A, 3*(ig-1)+2, 3*(ig-1)+2, LS_TGV);
		    }
	    }
    }
    
    for (k = 1; k <= mesh->nt; k++) {
	    ptt = &mesh->tria[k];
      ref =getRefb(mesh,ptt->ref);
	    if (ref) {
		    for (i = 0; i < 3; i++) {
		    	ig = ptt->v[i];
		        csrSet(A, 3*(ig-1)+0, 3*(ig-1)+0, LS_TGV);
		        csrSet(A, 3*(ig-1)+1, 3*(ig-1)+1, LS_TGV);
		        csrSet(A, 3*(ig-1)+2, 3*(ig-1)+2, LS_TGV);
		    }
	    }
    }

    return 1;
}

pCsr matA_P1_3d(pMesh mesh, pSol sol) {
    pCsr   A = 0;
    pTetra pt;
    double *a, *b, *c, *d, DeD[81], m[9], im[9], Ae[12][12], mm[9][12], nn[9][12], Dp[3][4];
    double lambda, mu, vol;
    int    i, j, k, s, ia, ja, il, ic, ig, jg, nr, nc, nbe;
    
    /* memory allocation (rough estimate) */
    nr  = nc = 3*mesh->np;
    nbe = 60*mesh->np;
    A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

    memset(DeD,0,81*sizeof(double));
    
    /* Dp */
    Dp[0][0]=1;  Dp[0][1]=0;  Dp[0][2]=0;  Dp[0][3]=-1;
    Dp[1][0]=0;  Dp[1][1]=1;  Dp[1][2]=0;  Dp[1][3]=-1;
    Dp[2][0]=0;  Dp[2][1]=0;  Dp[2][2]=1;  Dp[2][3]=-1;
    
    /* Fill stiffness matrix A */
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] )  continue;
        
        /* tD E D */
        if ( !getMat(sol,pt->ref,&lambda,&mu) )  continue;
        DeD[0]  = DeD[40] = DeD[80] = 2*mu + lambda;
        DeD[4]  = DeD[8]  = DeD[36] = DeD[44] = DeD[72] = DeD[76] = lambda;
        DeD[10] = DeD[12] = DeD[20] = DeD[24] = DeD[28] = DeD[30] = mu;
        DeD[50] = DeD[52] = DeD[56] = DeD[60] = DeD[68] = DeD[70] = mu;
        
        /* measure of K */
        a = &mesh->point[pt->v[0]].c[0];
        b = &mesh->point[pt->v[1]].c[0];
        c = &mesh->point[pt->v[2]].c[0];
        d = &mesh->point[pt->v[3]].c[0];
        
        /* mm = tB^-1 */
        for (i=0; i<3; i++) {
            m[i+0] = a[i] - d[i];
            m[i+3] = b[i] - d[i];
            m[i+6] = c[i] - d[i];
        }
        if ( !invmatg(m,im) )  return(0);
        
        vol = volume(a,b,c,d);
		assert(vol>0);
        
        /* mm = (tBt^-1) Dp */
        memset(mm,0,9*12*sizeof(double));
        for (i=0; i<3; i++) {
            for (j=0; j<4; j++) {
                for (s=0; s<3; s++)
                    mm[i][j]   += im[i*3+s] * Dp[s][j];
                mm[i+3][j+4] = mm[i][j];
                mm[i+6][j+8] = mm[i][j];
            }
        }
        
        /* nn = DeD mm */
        for (i=0; i<9; i++) {
            for (j=0; j<12; j++) {
                nn[i][j] = 0.0;
                for (s=0; s<9; s++)
                    nn[i][j] += DeD[i*9+s] * mm[s][j];
            }
        }
        
        /* Ae = vol tmm nn */
        memset(Ae,0,12*12*sizeof(double));
        for (i=0; i<12; i++) {
            for (j=i; j<12; j++) {
                for (s=0; s<9; s++)
                    Ae[i][j] += vol * mm[s][i] * nn[s][j];
            }
        }
        
        /* stifness matrix */
        for (i=0; i<12; i++) {
            ig = pt->v[i % 4];
            ia = 3*(ig-1) + (i / 4);
            for (j=i; j<12; j++) {
                if ( fabs(Ae[i][j]) < EPSD )  continue;
                jg = pt->v[j % 4];
                ja = 3*(jg-1) + (j / 4);
                if ( ia < ja ) {
                    il = ia;
                    ic = ja;
                }
                else {
                    il = ja;
                    ic = ia;
                }
                csrPut(A,il,ic,Ae[i][j]);
            }
        }
    }
    setTGV_3d(mesh,A);
    csrPack(A);

    return A;
}

/* find boundary conds in list */
pCl getCl(pSol sol, int ref, int elt) {
    pCl     pcl;
    int     i;
    
    for (i=0; i<sol->nbcl; i++) {
        pcl = &sol->cl[i];
        if ( (pcl->ref == ref) && (pcl->elt == elt) )  return(pcl);
    }
    return 0;
}

/* retrieve physical properties in list */
int getMat(pSol sol, int ref, double *lambda, double *mu) {
    pMat   pm;
    int    i;
    
    for (i=0; i<sol->nmat; i++) {
        pm = &sol->mat[i];
        if ( pm->ref == ref ) {
            *lambda = pm->lambda;
            *mu     = pm->mu;
            return(1);
        }
    }
  *lambda = LS_LAMBDA;
  *mu     = LS_MU;
  
 // *lambda = LS_E * LS_NU /( (1 + LS_NU)*(1-2*LS_NU));
 // *mu     = LS_E /( 2*(1 + LS_NU));

    return 1;
}


static int setTGV_2d(pMesh mesh, pCsr A) {
    pEdge    pa;
    pTria    ptt;
    int      k, ig, i, ref;

    for (k = 1; k <= mesh->nt; k++) {
	    ptt = &mesh->tria[k];
      ref =getRefel(mesh,ptt->ref);
      if (ref) {
        for (i = 0; i < 3; i++) {
		    	ig = ptt->v[i];
		        csrSet(A, 2*(ig-1)+0, 2*(ig-1)+0, LS_TGV);
		        csrSet(A, 2*(ig-1)+1, 2*(ig-1)+1, LS_TGV);
		    }
	    }
    }

    for (k=1; k<=mesh->na; k++) {
        pa = &mesh->edge[k];
      ref =getRefb(mesh,pa->ref);
       if (ref) {
        for(i=0;i<2;i++){
                ig = pa->v[i];
                csrSet(A,2*(ig-1)+0,2*(ig-1)+0,LS_TGV);
                csrSet(A,2*(ig-1)+1,2*(ig-1)+1,LS_TGV);
            }
        }
    }
    
    return 1;
}


pCsr matA_P1_2d(pMesh mesh, pSol sol) {
    pCsr     A;
    pTria    pt;
    double   Ae[6][6], DeD[4][4], m[2][2], mm[4][6], nn[4][6], *a, *b, *c;
    double   lambda, mu, det, idet, area;
    int      nr, nc, nbe, i, j, k, s, ia, ja, ig, jg, il, ic;
    
    fprintf(stdout,"\n  --  Begin building 2D matrix \n");
    
    /* memory allocation (rough estimate) */
    nr  = nc = 2*mesh->np;
    nbe = 40*mesh->np;
    A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);
    
    memset(mm,0,4*6*sizeof(double));
    memset(DeD,0,16*sizeof(double));
    
    /* store values in A */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] )  continue;
        
        /* tD E D */
        if ( !getMat(sol,pt->ref,&lambda,&mu) )  continue;
        DeD[0][0] = DeD[3][3] = 2*mu + lambda;
        DeD[0][3] = DeD[3][0] = lambda;
        DeD[1][1] = DeD[1][2] = DeD[2][1] = DeD[2][2] = mu;
        
        /* measure of K */
        a = &mesh->point[pt->v[0]].c[0];
        b = &mesh->point[pt->v[1]].c[0];
        c = &mesh->point[pt->v[2]].c[0];
        
        /* m = tBT^-1 */
        det  = (b[1]-c[1])*(a[0]-c[0])-(a[1]-c[1])*(b[0]-c[0]);
        if ( det < EPSD )  continue;
        idet = 1.0 / det;
        m[0][0] = idet*(b[1]-c[1]);    m[0][1] = idet*(c[1]-a[1]);
        m[1][0] = idet*(c[0]-b[0]);    m[1][1] = idet*(a[0]-c[0]);
        
        /* mm = (tBT^-1) Dp */
        mm[0][0] = mm[2][3] = m[0][0];
        mm[0][1] = mm[2][4] = m[0][1];
        mm[0][2] = mm[2][5] = -(m[0][0]+m[0][1]);
        mm[1][0] = mm[3][3] = m[1][0];
        mm[1][1] = mm[3][4] = m[1][1];
        mm[1][2] = mm[3][5] = -(m[1][0]+m[1][1]);
        
        /* nn = DeD mm */
        for (i=0; i<4; i++) {
            for (j=0; j<6; j++) {
                nn[i][j] = 0.0;
                for (s=0; s<4; s++)
                    nn[i][j] += DeD[i][s] * mm[s][j];
            }
        }
        
        area = 0.5 * det;
        assert(area>0);
        memset(Ae,0,6*6*sizeof(double));
        /* Ae = tmm * nn */
        for (i=0; i<6; i++) {
            for (j=i; j<6; j++) {
                for (s=0; s<4; s++)
                    Ae[i][j] += area * mm[s][i] * nn[s][j];
            }
        }
        
        /* stifness matrix */
        for (i=0; i<6; i++) {
            ig = pt->v[i % 3];
            ia = 2*(ig-1) + (i / 3);
            for (j=i; j<6; j++) {
                if ( fabs(Ae[i][j]) < EPSD )  continue;
                jg = pt->v[j % 3];
                ja = 2*(jg-1) + (j / 3);
                if ( ia < ja ) {
                    il = ia;
                    ic = ja;
                }
                else {
                    il = ja;
                    ic = ia;
                }
                csrPut(A,il,ic,Ae[i][j]);
            }
        }
    }
    
    setTGV_2d(mesh,A);
 
    csrPack(A);
    //if ( abs(info.imprim) > 5 || info.ddebug )
       // fprintf(stdout,"     A: %6d x %6d  sparsity %7.4f%%\n",nr,nc,100.0*A->nbe/nr/nc);
    return A;
}

/* linear elasticity */
int elas(pInstance instance, double *F) {

    double   err;
    int      ier, nit;
    pCsr 	 A = 0;
  
    err = LS_RES;
    nit = 1000;
    A = matA_P1(&instance->mesh_omega0,&instance->sol_omega0);
	assert(A != 0);

    ier = csrPrecondGrad(A,instance->sol_omega0.u,F,&err,&nit,0);
    if  ( info.ncpu > 1 )  csrStop();
    if ( ier <= 0 )  fprintf(stdout,"          ## SOL NOT CONVERGED: ier = %d  err = %E  nit = %d \n",ier,err,nit);
    csrFree(A);
    return (ier > 0);
}
