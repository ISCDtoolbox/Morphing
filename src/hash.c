#include "morphing.h"

#define KA     31
#define KB     57
#define KC     79
#define KTA     7
#define KTB    11

unsigned char idir[5]     = {0,1,2,0,1};
unsigned char idirt[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };

int hashelt_3d(pMesh mesh) {
    pTetra    pt,pt1;
    int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
    int      *hcode,*link,inival,hsize;
    unsigned char   i,ii,i1,i2,i3;
    unsigned int    key;
    
    /* memory alloc */
    hcode = (int*)calloc(mesh->ne+1,sizeof(int));
    assert(hcode);
    link  = mesh->adja;
    hsize = mesh->ne;
    
    /* init */
    inival = 2147483647;
    for (k=0; k<=mesh->ne; k++)
        hcode[k] = -inival;
    
    /* build hash table */
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] )  continue;
        for (i=0; i<4; i++) {
            i1 = idirt[i][0];
            i2 = idirt[i][1];
            i3 = idirt[i][2];
            mins = LS_MIN(pt->v[i1],pt->v[i2]);
            mins = LS_MIN(mins,pt->v[i3]);
            maxs = LS_MAX(pt->v[i1],pt->v[i2]);
            maxs = LS_MAX(maxs,pt->v[i3]);
            
            /* compute key */
            sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
            key = KA*mins + KB*maxs + KC*sum;
            key = key % hsize + 1;
            /* insert */
            iadr = 4*(k-1) + i +1 ;
            link[iadr] = hcode[key];
            hcode[key] = -iadr;
        }
    }
    /* set adjacency */
    for (l=4*mesh->ne; l>0; l--) {
        if ( link[l] >= 0 )  continue;
        k = ((l-1) >> 2) + 1;
        i = (l-1) % 4;
        i1 = idirt[i][0];
        i2 = idirt[i][1];
        i3 = idirt[i][2];
        pt = &mesh->tetra[k];
        
        sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];
        mins = LS_MIN(pt->v[i1],pt->v[i2]);
        mins = LS_MIN(mins,pt->v[i3]);
        maxs = LS_MAX(pt->v[i1],pt->v[i2]);
        maxs = LS_MAX(maxs,pt->v[i3]);
        
        /* accross link */
        ll = -link[l];
        pp = 0;
        link[l] = 0;
        while ( ll != inival ) {
            kk = ((ll-1) >> 2) + 1;
            ii = (ll-1) % 4;
            i1 = idirt[ii][0];
            i2 = idirt[ii][1];
            i3 = idirt[ii][2];
            pt1  = &mesh->tetra[kk];
            sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
            if ( sum1 == sum ) {
                mins1 = LS_MIN(pt1->v[i1],pt1->v[i2]);
                mins1 = LS_MIN(mins1,pt1->v[i3]);
                if ( mins1 == mins ) {
                    maxs1 = LS_MAX(pt1->v[i1],pt1->v[i2]);
                    maxs1 = LS_MAX(maxs1,pt1->v[i3]);
                    if ( maxs1 == maxs ) {
                        /* adjacent found */
                        if ( pp != 0 )  link[pp] = link[ll];
                        link[l] = 4*kk + ii;
                        link[ll]= 4*k + i;
                        break;
                    }
                }
            }
            pp = ll;
            ll = -link[ll];
        }
    }
    
    free(hcode);
    return(1);
}


int hashelt_2d(pMesh mesh) {
    pTria     pt,pt1;
    int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,iadr;
    int      *hcode,*link,inival,hsize;
    unsigned char  *hvoy,i,ii,i1,i2;
    unsigned int    key;
    
    /* memory alloc */
    hcode = (int*)calloc(mesh->nt+1,sizeof(int));
    assert(hcode);
    link  = mesh->adja;
    hsize = mesh->nt;
    hvoy  = (unsigned char*)hcode;
    
    /* init */
    inival = 2147483647;
    for (k=0; k<=mesh->nt; k++)
        hcode[k] = -inival;
    
    /* build hash table */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] )  continue;
        
        for (i=0; i<3; i++) {
            i1 = idir[i+1];
            i2 = idir[i+2];
            if ( pt->v[i1] < pt->v[i2] ) {
                mins = pt->v[i1];
                maxs = pt->v[i2];
            }
            else {
                mins = pt->v[i2];
                maxs = pt->v[i1];
            }
            
            /* compute key */
            key = KTA*mins + KTB*maxs;
            key = key % hsize + 1;
            
            /* insert */
            iadr = 3*(k-1) + i+1;
            link[iadr] = hcode[key];
            hcode[key] = -iadr;
        }
    }
    
    /* set adjacency */
    for (l=3*mesh->nt; l>0; l--) {
        if ( link[l] >= 0 )  continue;
        k = (l-1) / 3 + 1;
        i = (l-1) % 3;
        i1 = idir[i+1];
        i2 = idir[i+2];
        pt = &mesh->tria[k];
        
        mins = LS_MIN(pt->v[i1],pt->v[i2]);
        maxs = LS_MAX(pt->v[i1],pt->v[i2]);
        
        /* accross link */
        ll = -link[l];
        pp = 0;
        link[l] = 0;
        hvoy[l] = 0;
        while ( ll != inival ) {
            kk = (ll-1) / 3 + 1;
            ii = (ll-1) % 3;
            i1 = idir[ii+1];
            i2 = idir[ii+2];
            pt1  = &mesh->tria[kk];
            if ( pt1->v[i1] < pt1->v[i2] ) {
                mins1 = pt1->v[i1];
                maxs1 = pt1->v[i2];
            }
            else {
                mins1 = pt1->v[i2];
                maxs1 = pt1->v[i1];
            }
            
            if ( mins1 == mins  && maxs1 == maxs ) {
                /* adjacent found */
                if ( pp != 0 )  link[pp] = link[ll];
                link[l] = 3*kk + ii;
                link[ll]= 3*k + i;
                break;
            }
            pp = ll;
            ll = -link[ll];
        }
    }
    
    free(hcode);
    return(1);
}

