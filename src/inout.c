#include "morphing.h"

extern Info   info;

/* read mesh */
int loadMesh(pMesh mesh) {
    pPoint       ppt;
    pTria        pt1;
    pTetra       ptt;
    pEdge       pa;
    float        fp1,fp2,fp3;
    int          k,inm,i;
    char        *ptr,data[256];
  
    mesh->dim=0;
  
    strcpy(data,mesh->name);
    ptr = strstr(data,".mesh");
    if ( !ptr ) {
        strcat(data,".mesh");
        if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
            strcat(data,"b");
            if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
                fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
                return(0);
            }
        }
    }
    else {
    	if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
           fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
     	   return(0);
    	}
    }
    fprintf(stdout,"  %%%% %s OPENED\n",data);
    
    if ( abs(info.imprim) > 3 )
        fprintf(stdout,"  -- READING DATA FILE %s\n",data);
    
    mesh->np = GmfStatKwd(inm,GmfVertices);
    mesh->nt = GmfStatKwd(inm,GmfTriangles);
    mesh->ne = GmfStatKwd(inm,GmfTetrahedra);
    mesh->na = GmfStatKwd(inm,GmfEdges);
    
    if ( !mesh->np ) {
        fprintf(stdout,"  ** MISSING DATA\n");
        return(0);
    }
    
    /* memory alloc */
    mesh->point = (pPoint)calloc(mesh->np+1,sizeof(Point));
    assert(mesh->point);
    if ( mesh->ne ) {
        mesh->tetra  = (pTetra)calloc(mesh->ne+1,sizeof(Tetra));
        assert(mesh->tetra);
    }
    if ( mesh->nt ) {
        mesh->tria  = (pTria)calloc(mesh->nt+1,sizeof(Tria));
        assert(mesh->tria);
    }
    
    if ( mesh->dim == 2 ) {
        
        if ( mesh->na ) {
            mesh->edge  = (pEdge)calloc(mesh->na+1,sizeof(Edge));
            assert(mesh->edge);
        }
        /* read mesh vertices */
        GmfGotoKwd(inm,GmfVertices);
        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            if ( mesh->ver == GmfFloat ) {
                GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ppt->ref);
                ppt->c[0] = fp1;
                ppt->c[1] = fp2;
            }
            else
                GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->ref);
        }
        /* read mesh triangles */
        GmfGotoKwd(inm,GmfTriangles);
        for (k=1; k<=mesh->nt; k++) {
            pt1 = &mesh->tria[k];
            GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
            for (i=0; i<3; i++) {
                ppt = &mesh->point[pt1->v[i]];
                pt1->g[0] += ppt->c[0];
                pt1->g[1] += ppt->c[1];
                if ( !ppt->s )  ppt->s = k;
            }
            pt1->g[0] /= 3.0;
            pt1->g[1] /= 3.0;
        }
        
        /* read mesh edges */
        GmfGotoKwd(inm,GmfEdges);
        for (k=1; k<=mesh->na; k++) {
            pa = &mesh->edge[k];
            if ( mesh->ver == GmfFloat ) {
                GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
                
            }
            else
                GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
        }

        if ( abs(info.imprim) > 4 ) {
            fprintf(stdout,"  %%%% NUMBER OF VERTICES  %8d\n",mesh->np);
            if ( mesh->na )  fprintf(stdout,"  %%%% NUMBER OF EDGES     %8d\n",mesh->na);
            if (mesh->nt)    fprintf(stdout,"  %%%% NUMBER OF TRIANGLES %8d\n",mesh->nt);
        }
    }
    
    else {
        GmfGotoKwd(inm,GmfVertices);
        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            if ( mesh->ver == GmfFloat ) {
                GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
                ppt->c[0] = fp1;
                ppt->c[1] = fp2;
                ppt->c[2] = fp3;
            }
            else
                GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
        }
        /* read triangles */
        GmfGotoKwd(inm,GmfTriangles);
        for (k=1; k<=mesh->nt; k++) {
            pt1 = &mesh->tria[k];
            GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
        }
        /* read tetrahedra */
        GmfGotoKwd(inm,GmfTetrahedra);
        for (k=1; k<=mesh->ne; k++) {
            ptt = &mesh->tetra[k];
            GmfGetLin(inm,GmfTetrahedra,&ptt->v[0],&ptt->v[1],&ptt->v[2],&ptt->v[3],&ptt->ref);
            for (i=0; i<4; i++) {
                ppt = &mesh->point[ptt->v[i]];
                ptt->g[0] += ppt->c[0];
                ptt->g[1] += ppt->c[1];
                ptt->g[2] += ppt->c[2];
                if ( !ppt->s )  ppt->s = k;
            }
            ptt->g[0] *= 0.25;
            ptt->g[1] *= 0.25;
            ptt->g[2] *= 0.25;
        }
        
        if ( abs(info.imprim) > 4 ) {
            fprintf(stdout,"  %%%% NUMBER OF VERTICES   %8d\n",mesh->np);
            if ( mesh->nt )  fprintf(stdout,"  %%%% NUMBER OF TRIANGLES  %8d\n",mesh->nt);
            if ( mesh->ne )  fprintf(stdout,"  %%%% NUMBER OF TETRAHEDRA %8d\n",mesh->ne);
        }
    }
    
    GmfCloseMesh(inm);
    return(1);
}

/* write mesh */
int saveMesh(pMesh mesh, pMesh mesh2, int n) {
    pPoint       ppt;
    pTria        pt1;
    pEdge        pa;
    pTetra       pt;
    int          k,inm;
    char         *ptr,data[128];
    
  mesh->ver = GmfDouble;
  strcpy(data,mesh2->name);
  ptr = strstr(data,".mesh");
  if ( ptr ) *ptr = '\0';
  sprintf(data,"%s.%d.mesh",data,n);
  
  if ( !(inm = GmfOpenMesh(data,GmfWrite,mesh->ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
    fprintf(stdout,"  %%%% %s OPENED\n",data);
    
    GmfSetKwd(inm,GmfVertices,mesh->np);
    if ( mesh->dim == 2 ) {
        for(k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            GmfSetLin(inm,GmfVertices,ppt->c[0],ppt->c[1],ppt->ref);
        }
    }
    else {
        for(k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            
            GmfSetLin(inm,GmfVertices,ppt->c[0],ppt->c[1],
                      ppt->c[2],ppt->ref);
        }
    }
    /* write triangles */
    GmfSetKwd(inm,GmfTriangles,mesh->nt);
    for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        GmfSetLin(inm,GmfTriangles,pt1->v[0],pt1->v[1],pt1->v[2],pt1->ref);
    }
    
    /* write tetrahedra */
    GmfSetKwd(inm,GmfTetrahedra,mesh->ne);
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] )  continue;
        GmfSetLin(inm,GmfTetrahedra,pt->v[0],pt->v[1],pt->v[2],pt->v[3],pt->ref);
    }
  

  /* write edges */
  if(mesh->dim==2){
  GmfSetKwd(inm,GmfEdges,mesh->na);
  for (k=1; k<=mesh->na; k++) {
    pa  = &mesh->edge[k];
    if ( !pa->v[0] )  continue;
    GmfSetLin(inm,GmfEdges,pa->v[0],pa->v[1],pa->ref);
  }
  }
  
    GmfCloseMesh(inm);
    return(1);
}

/* load scalar solution field */
int loadSol(pSol sol) {
    float       buf[GmfMaxTyp];
    double      bufd[GmfMaxTyp];
    int         i,k,type,inm,typtab[GmfMaxTyp],offset;
    char       *ptr,data[128];
 
    if ( !sol->name )  return(-1);
    strcpy(data,sol->name);
    ptr = strstr(data,".mesh");

 	if ( ptr ) *ptr = '\0';
 	strcat(data,".sol");
        if ( !(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
            strcat(data,"b");
            if ( !(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
                fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
                return(0);
            }
        }
    fprintf(stdout,"  %%%% %s OPENED\n",data);
    
    if ( abs(info.imprim) > 3 )
        fprintf(stdout,"  -- READING DATA FILE %s\n",data);
    
    sol->np = GmfStatKwd(inm,GmfSolAtVertices,&type,&offset,&typtab);
    fprintf(stdout,"  -- READING DATA FILE %d %d %d \n",type,offset,typtab[0]);
    if ( !sol->np  )  return(-1);
    sol->size[0] = offset;
    
    if (typtab[0]==2) {
        
        /* alloc */
        sol->u  = (double*)calloc(sol->dim*sol->np,sizeof(double));
        assert(sol->u);
        
        /* read mesh solutions */
        GmfGotoKwd(inm,GmfSolAtVertices);
        if ( sol->ver == GmfFloat ) {
            for (k=0; k<sol->np; k++) {
                GmfGetLin(inm,GmfSolAtVertices,&buf);
                for (i=0; i<sol->dim; i++)
                    sol->u[sol->dim*k+i] = buf[i];
            }
        }
        else {
            for (k=0; k<sol->np; k++) {
                GmfGetLin(inm,GmfSolAtVertices,bufd);
                for (i=0; i<sol->dim; i++)
                    sol->u[sol->dim*k+i] = bufd[i];
            }
        }
    }
    else if (typtab[0]==1) {
        sol->u  = (double*)calloc(sol->np,sizeof(double));
        assert(sol->u);
        sol->valp1  = (double*)calloc(sol->np,sizeof(double));
        assert(sol->valp1 );
        
        GmfGotoKwd(inm,GmfSolAtVertices);
        if ( sol->ver == GmfFloat ) {
            for (k=0; k<sol->np; k++) {
                GmfGetLin(inm,GmfSolAtVertices,&buf);
                //sol->u[k] = buf[0];
                sol->valp1[k] = buf[0];
            }
        }
        else {
            for (k=0; k<sol->np; k++) {
                GmfGetLin(inm,GmfSolAtVertices,bufd);
                // sol->u[k] = bufd[0];
                sol->valp1[k] = bufd[0];
            }
        }
    }
 
    GmfCloseMesh(inm);
    return(1);
}

/* save vectorial solution field */
int saveSol(pSol sol, pMesh mesh, pMesh mesh2, int n) {
    double       dbuf[GmfMaxTyp];
    float        fbuf[GmfMaxTyp];
    int          k,ia,i,inm,type,typtab[GmfMaxTyp];
    char        *ptr,data[128];
  
  	strcpy(data,mesh2->name);
  	ptr = strstr(data,".mesh");

 	 if ( ptr ) *ptr = '\0';
  	 sprintf(data,"%s.%d.depl.sol",data,n);

    if ( !(inm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
        return(0);
    }
    fprintf(stdout,"  %%%% %s OPENED\n",data);
    
    type = 1;
    typtab[0] = GmfVec;
    
    /* write sol */
    GmfSetKwd(inm,GmfSolAtVertices,sol->np,type,typtab);
    if ( sol->ver == GmfFloat ) {
        for (k=0; k<sol->np; k++) {
            ia = sol->dim*k;
            for (i=0; i<sol->dim; i++)
                fbuf[i] = sol->p[ia+i];
            GmfSetLin(inm,GmfSolAtVertices,fbuf);
        }
    }
    else {
        for (k=0; k<sol->np; k++) {
            ia = sol->dim*k;
            for (i=0; i<sol->dim; i++)
                dbuf[i] = sol->p[ia+i];
            GmfSetLin(inm,GmfSolAtVertices,dbuf);
        }
    }
    
    GmfCloseMesh(inm);
    return(1);
}

/* save scalar signed distance function */
int saveDistance(pSol sol, pMesh mesh, pMesh mesh2, int n) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,i,inm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128];
  
  
  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  
   if ( ptr ) *ptr = '\0';
   sprintf(data,"%s.%d.sol",data,n);


  if ( !(inm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  
  type = 1;
  typtab[0] = GmfSca;
  
  /* write sol */
  GmfSetKwd(inm,GmfSolAtVertices,sol->np,type,typtab);
  if ( sol->ver == GmfFloat ) {
    for (k=0; k<sol->np; k++) {
      for (i=0; i<typtab[0]; i++)
      fbuf[i] = sol->d[k];
      GmfSetLin(inm,GmfSolAtVertices,fbuf);
   }
  }
  else {
    for (k=0; k<sol->np; k++) {
      for (i=0; i<typtab[0]; i++)
        dbuf[i] = sol->d[k];
      GmfSetLin(inm,GmfSolAtVertices,dbuf);
      }
    }
  
  GmfCloseMesh(inm);
  return(1);
}



/* superopose the contours of mesh and mesh2 */
int saveContour(pMesh mesh, pMesh mesh2) {
  pPoint       ppt;
  pTria        pt1;
  pEdge        pa;
  int          k,inm;
  char        *ptr,data[128];
  
  mesh->ver = GmfDouble;
  
  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".fin.mesh");
    if( !(inm = GmfOpenMesh(data, GmfWrite, mesh->ver,mesh->dim)) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        return(0);
    }
  }
  else {
  	*ptr = '\0';
  	if( !(inm = GmfOpenMesh(data, GmfWrite, mesh->ver,mesh->dim)) ) {
    	fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
    	return(0);
  	}
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  
  GmfSetKwd(inm,GmfVertices,mesh->np+ mesh2->np);
  if ( mesh->dim == 2 ) {
    for(k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      GmfSetLin(inm,GmfVertices,ppt->c[0],ppt->c[1],ppt->ref);
    }
    for(k=1; k<=mesh2->np; k++) {
      ppt = &mesh2->point[k];
      GmfSetLin(inm,GmfVertices,ppt->c[0],ppt->c[1],ppt->ref);
    }
    
  }
  else {
    for(k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      
      GmfSetLin(inm,GmfVertices,ppt->c[0],ppt->c[1],
                ppt->c[2],ppt->ref);
    }
  }
  
  /* write triangles */
  GmfSetKwd(inm,GmfTriangles,mesh2->nt);
  for (k=1; k<=mesh2->nt; k++) {
    pt1 = &mesh2->tria[k];
    GmfSetLin(inm,GmfTriangles,pt1->v[0]+ mesh2->np,pt1->v[1]+mesh2->np,pt1->v[2]+mesh2->np,pt1->ref);
  }
  
  /* write edges */
  GmfSetKwd(inm,GmfEdges,mesh->na + mesh2->na);
  for (k=1; k<=mesh->na; k++) {
    pa  = &mesh->edge[k];
    if ( pa->ref == 10 ) {
      GmfSetLin(inm,GmfEdges,pa->v[0],pa->v[1],7);
    }
  }
  
  for (k=1; k<=mesh2->na; k++) {
    pa  = &mesh2->edge[k];
    if ( pa->ref == 1) {
      GmfSetLin(inm,GmfEdges,pa->v[0]+mesh->np,pa->v[1]+mesh->np,1);
    }
  }
  
  GmfCloseMesh(inm);
  return(1);
}




