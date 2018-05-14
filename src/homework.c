
#include"fem.h"
/*
Louis NAVARRE : 1235-16-00
Nahi NASSAR : 1269-16-00

Projet FEM 2018
*/
#define VEXT 7.0


# ifndef NOCOUETTECREATE

femCouetteProblem *femCouetteCreate(const char *filename, femGrains *grains)
{
    femCouetteProblem *theProblem = malloc(sizeof(femCouetteProblem));
    theProblem->mesh              = femMeshRead(filename);           
    theProblem->edges             = femEdgesCreate(theProblem->mesh);  
    theProblem->space             = femDiscreteCreate(3,FEM_TRIANGLE);
    theProblem->rule              = femIntegrationCreate(3,FEM_TRIANGLE); 
    theProblem->system            = femFullSystemCreate(theProblem->mesh->nNode);
    theProblem->system2           = femFullSystemCreate(theProblem->mesh->nNode);
    theProblem->mu                = 1.0;
    theProblem->grains            = grains;
    theProblem->norm              = malloc(sizeof(double)*theProblem->system->size);
    return theProblem;
}

# endif
# ifndef NOCOUETTEFREE

void femCouetteFree(femCouetteProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femFullSystemFree(theProblem->system2);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    femGrainsFree(theProblem->grains);
    free(theProblem);
}
    

# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
  //node est le noeud du milieu
  int node = theMesh->nLocalNode;
  for(int j=0;j<node;j++)
  {
    //on remplit map, X et Y pour node
    map[j] = theMesh->elem[i*node+j];
    x[j]   = theMesh->X[map[j]];
    y[j]   = theMesh->Y[map[j]];
  }
}

# endif

void femIterativeSolverAssemble(femIterativeSolver* mySolver, femFullSystem *theSystem)
{
    int i,j,myRow;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    double *s  = mySolver->S;
    double *d  = mySolver->D;
    double *r  = mySolver->R;
    int iter   = mySolver->iter;
    /* 
    * On ne calcule R qu'a la premiere iteration, car r(k+1) est
    * calcule a chaque iteration
    * S, lui, doit etre calcule a chaque iteration car la direction
    * d change a chaque iteration
    */
    for (i = 0; i < mySolver->size; ++i)
    {
      for (j = 0; j < mySolver->size; ++j)
      {
        if (iter == 0)
        {
          r[i] += A[i][j] * mySolver->X[j];
        }
        s[i] += A[i][j] * d[j];
      }
      if (iter == 0)
      {
        r[i] -= B[i];
      }
    }
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) 
{   
    mySolver->R[myNode] = myValue;
    mySolver->S[myNode] = myValue;
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i]) * (mySolver->R[i]);
    }

    int taille = mySolver->size;
    double *s  = mySolver->S;
    double *d  = mySolver->D;
    double *r  = mySolver->R;
    double *dx = mySolver->X;

    if(mySolver->iter == 1){
        for(i=0;i<taille;i++){
            d[i]  = r[i];
            dx[i] = 0;
        }
    }
    else{
        double da = 0.0;
        for(i = 0; i < taille; i++){
            da += s[i]*r[i];
        }

        double alpha = -error / da;
        double nb = 0.0;
        for(i = 0; i < taille; i++){
            dx[i] += alpha*d[i];
            r[i]   = r[i] + alpha * s[i];
            nb    += r[i] * r[i];
        }

        double beta = nb / error;
        for(i = 0; i < taille; i++){
            d[i] = r[i] + beta * d[i];
            s[i] = 0.0;
        }
    }
    mySolver->error = sqrt(error);

    return mySolver->X;
}

# ifndef FEMCOUETTE

void femCouetteSolve(femCouetteProblem *theProblem){
  femIterativeSolver *theSolver  = femIterativeSolverCreate(theProblem->system->size);
  femIterativeSolver *theSolver2 = femIterativeSolverCreate(theProblem->system2->size);

  int testconvergence = 0, testconvergence2 = 0;

  do{
    femIterativeSolverAssemble(theSolver, theProblem->system);
    femIterativeSolverEliminate(theSolver);

    femIterativeSolverAssemble(theSolver2, theProblem->system2);
    femIterativeSolverEliminate(theSolver2);

    testconvergence  = femIterativeSolverConverged(theSolver);
    testconvergence2 = femIterativeSolverConverged(theSolver2);
  }while(testconvergence == 0 && testconvergence2 == 0);

  if(testconvergence == -1 || testconvergence == -1){
    printf("Trop d'iterations\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0; i<theSolver->size;i++){
    theProblem->system->B[i]  = theSolver->X[i];
    theProblem->system2->B[i] = theSolver2->X[i];
    theProblem->norm[i]       = sqrt(theProblem->system->B[i]*theProblem->system->B[i] + theProblem->system2->B[i]*theProblem->system2->B[i]);
  }

  free(theSolver);
  free(theSolver2);
}

#endif

void isomorphisme(double x, double y, double *x_loc, double *y_loc, double *to_return){
  double x_iso, y_iso;

  x_iso = -(x*y_loc[0] - x_loc[0]*y - x*y_loc[2] + x_loc[2]*y + x_loc[0]*y_loc[2] - x_loc[2]*y_loc[0])/(x_loc[0]*y_loc[1] - x_loc[1]*y_loc[0] - x_loc[0]*y_loc[2] + x_loc[2]*y_loc[0] + x_loc[1]*y_loc[2] - x_loc[2]*y_loc[1]);
  y_iso =  (x*y_loc[0] - x_loc[0]*y - x*y_loc[1] + x_loc[1]*y + x_loc[0]*y_loc[1] - x_loc[1]*y_loc[0])/(x_loc[0]*y_loc[1] - x_loc[1]*y_loc[0] - x_loc[0]*y_loc[2] + x_loc[2]*y_loc[0] + x_loc[1]*y_loc[2] - x_loc[2]*y_loc[1]);

  to_return[0] = x_iso;
  to_return[1] = y_iso;
}

# ifndef FINDELEMENT

void findElement(femCouetteProblem *theProblem){
  femGrains *theGrains = theProblem->grains;
  femMesh *theMesh     = theProblem->mesh;
  int n_g              = theGrains->n;
  int k;
  int *elem_g          = theGrains->elem;
  double *x_g          = theGrains->x;
  double *y_g          = theGrains->y;

  double *X_m          = theMesh->X;
  double *Y_m          = theMesh->Y;
  int *elem_m          = theMesh->elem;
  int n_e              = theMesh->nElem;
  double *B            = theProblem->system->B;
  double *B2           = theProblem->system2->B;

  int i,j;
  double x_elem[3], y_elem[3];
  int map[3];
  double iso[2];

  for(i = 0; i < n_g; i++){
    double x_loc = x_g[i];
    double y_loc = y_g[i];
    int dedans = 0;

    for(j = 0; j < n_e && dedans == 0; j++){
      femMeshLocal(theMesh, j, map, x_elem, y_elem);
      isomorphisme(x_loc, y_loc, x_elem, y_elem, iso);

      if(iso[0] <= 1 && iso[0] >= 0 && iso[1] <= 1 && iso[1] >= 0 && iso[1] + iso[0] <= 1){
        dedans = 1;
        elem_g[i] = j;
      }
    }

  }

}

# endif

# ifndef NOPOISSONASSEMBLE


void femCouetteAssemble(femCouetteProblem *theProblem)
{
  femGrains *theGrains = theProblem->grains;
  findElement(theProblem);
  int option = 1;
  femMesh *theMesh          = theProblem->mesh;
  femEdges *theEdges        = theProblem->edges;
  femDiscrete *theSpace     = theProblem->space;
  femIntegration *theRule   = theProblem->rule;
  femFullSystem *theSystem  = theProblem->system;
  femFullSystem *theSystem2 = theProblem->system2;
  double gamma              = theGrains->gamma;
  double radiusOut          = theGrains->radiusOut;
  double mu                 = theProblem->mu;

  double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
  int i,j,k,l;
  int map[4];

  for (i = 0; i < theMesh->nElem; i++)
  {
    //creation d'un maillage local
    femMeshLocal(theMesh,i,map,x,y);

    int nSpace=theSpace->n;

    for (j=0; j < 3; j++)
    {
      //valeurs d'integration
      double xsi    = theRule->xsi[j];
      double eta    = theRule->eta[j];
      double weight = theRule->weight[j];

      //fonctions de forme locales
      femDiscretePhi2(theSpace,xsi,eta,phi);
      //derivees fonctions de forme locales par rapport à xsi et eta
      femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

      //passage du plan xsi-eta au plan x-y
      double dxdxsi = 0, dxdeta = 0, dydxsi = 0, dydeta = 0 ;
      for (k = 0; k < 3; k++)
      {
        dxdxsi = dxdxsi + x[k]*dphidxsi[k];  
        dxdeta = dxdeta + x[k]*dphideta[k];   
        dydxsi = dydxsi + y[k]*dphidxsi[k];   
        dydeta = dydeta + y[k]*dphideta[k];
      }

      //calcul du jacobien local (voir slide 17 CM4)
      double J_e = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

      for (k = 0; k < nSpace; k++)
      {
        dphidx[k] = (dphidxsi[k] * dydeta - dphideta[k] * dydxsi) / J_e;
        dphidy[k] = (dphideta[k] * dxdxsi - dphidxsi[k] * dxdeta) / J_e;
      }

      for (k = 0; k < 3; k++)
      {
        for(l = 0; l < 3; l++)
        {
          double numIntegrateA          = (dphidx[k] * dphidx[l] + dphidy[k] * dphidy[l]) * J_e * weight;
          theSystem->A[map[k]][map[l]]  = theSystem->A[map[k]][map[l]] + mu*numIntegrateA;
          theSystem2->A[map[k]][map[l]] = theSystem2->A[map[k]][map[l]] + mu*numIntegrateA;
        }
        theSystem->B[map[k]]  = 0;
        theSystem2->B[map[k]] = 0;
      }
    }
  }

  int *elem = theGrains->elem;
  double xx[2];
  double iso[2];
  double myPhi[3];
  for(i = 0; i < theGrains->n; i++){
    int myElem = elem[i];
    if(myElem > 0){
      xx[0] = theGrains->x[i];
      xx[1] = theGrains->y[i];
      femMeshLocal(theMesh, myElem, map, x, y);
      isomorphisme(xx[0], xx[1], x, y, iso);
      _phi(iso[0], iso[1], myPhi);
      for(k = 0; k < 3; k++){
        for(l = 0; l < 3; l++){
          double grainsA                 = myPhi[k]*myPhi[l];
          theSystem->A[map[k]][map[l]]  += gamma*grainsA;
          theSystem2->A[map[k]][map[l]] += gamma*grainsA;
        }
        theSystem->B[map[k]]  = gamma*myPhi[k]*theGrains->vx[i];
        theSystem2->B[map[k]] = gamma*myPhi[k]*theGrains->vy[i];
      }
    }


  }

  int done = 0;
  for(i = 0; i < theEdges->nEdge && !done; i++){
    if(theEdges->edges[i].elem[1] == -1){
        double xe1 = theMesh->X[theEdges->edges[i].node[0]];
        double ye1 = theMesh->Y[theEdges->edges[i].node[0]];
        double xe2 = theMesh->X[theEdges->edges[i].node[1]];
        double ye2 = theMesh->Y[theEdges->edges[i].node[1]];
        double r1  = sqrt(xe1*xe1 + ye1*ye1);
        if(r1 < radiusOut*0.95){
          femFullSystemConstrain(theSystem,  theEdges->edges[i].node[0], 0); 
          femFullSystemConstrain(theSystem2, theEdges->edges[i].node[0], 0);
          femFullSystemConstrain(theSystem,  theEdges->edges[i].node[1], 0); 
          femFullSystemConstrain(theSystem2, theEdges->edges[i].node[1], 0);
        }
        else{
          double t1 = atan2(xe1,ye1); double t2 = atan2(xe2,ye2);
          femFullSystemConstrain(theSystem,  theEdges->edges[i].node[0],  VEXT*cos(t1)); 
          femFullSystemConstrain(theSystem2, theEdges->edges[i].node[0], -VEXT*sin(t1));
          femFullSystemConstrain(theSystem,  theEdges->edges[i].node[1],  VEXT*cos(t2)); 
          femFullSystemConstrain(theSystem2, theEdges->edges[i].node[1], -VEXT*sin(t2));
        }
    }
  }

  femCouetteSolve(theProblem);
}

# endif

# ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  
{
    int n = myGrains->n; 
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *r          = myGrains->r;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn         = myGrains->radiusIn;
    double rOut        = myGrains->radiusOut;
    
    double zeta = 0.0;

    int i,j,iContact = 0;

    /* Reinitialisation des corrections de vitesse */
    double dv;


    /* Correction de vitesse */

    /* Calcul d une correction de vitesse pour tenir compte des collisions entre les grains */
    for(i=0;i<n-1;++i){
        for(j=i+1;j<n;++j){
            double nx = x[j]-x[i];
            double ny = y[j]-y[i];
            double nn = sqrt(nx*nx+ny*ny);
            nx = nx/nn;
            ny = ny/nn;     
            if(iter==0){
                dv = dvContacts[iContact];
            }
            else{       
            double gamma          = nn-r[i]-r[j];
            double vn             =(vx[i]*nx+vy[i]*ny)-(vx[j]*nx+vy[j]*ny);
            dv                    =fmax(0,vn+dvContacts[iContact]-gamma/dt)-dvContacts[iContact];
            dvContacts[iContact] += dv;
            }

            vx[i] += -dv*nx*m[j]/(m[i]+m[j]);
            vy[i] += -dv*ny*m[j]/(m[i]+m[j]);
            vx[j] += dv*nx*m[i]/(m[i]+m[j]);
            vy[j] += dv*ny*m[i]/(m[i]+m[j]);

            zeta=fmax(zeta,fabs(dv));

            iContact++;
        }
    }
    
    /* Calcul d une correction de vitesse pour tenir compte des collisions avec la frontiere */
    for(i=0;i<n;i++){

        double norm = sqrt((x[i]*x[i])+(y[i]*y[i]));
        double nx = x[i]/norm;
        double ny = y[i]/norm;
        if(iter==0){
            dv = dvBoundary[i];
        }
        else {
        double vn         = (vx[i]*nx+vy[i]*ny);
        double gammaOut   = rOut-norm-r[i];
        double gammaIn    = norm-rIn-r[i];
        double dv1        = fmax(0,vn+dvBoundary[i]-gammaOut/dt);
        double dv2        = fmax(0,-vn-dvBoundary[i]-gammaIn/dt);
        dv                = dv1 - dv2-dvBoundary[i];
        dvBoundary[i]    += dv;
        }

        vx[i] += -dv*nx;
        vy[i] += -dv*ny;

        zeta=fmax(zeta,fabs(dv));

    }
    return zeta;
}

# endif

# ifndef REINITI

void reinitialiseMatrice(femCouetteProblem *theProblem){
  int i,j;
  for(i=0;i<theProblem->mesh->nNode;i++){
    for(j=0;j<theProblem->mesh->nNode;j++){
      theProblem->system->A[i][j]   = 0;
      theProblem->system2->A[i][j]  = 0;
    }
    theProblem->system->B[i]  = 0;
    theProblem->system2->B[i] = 0;
  }
}

#endif

# ifndef NOUPDATE


void femGrainsUpdate(femCouetteProblem *theProblem, double dt, double tol, double iterMax)
{
    femGrains *myGrains = theProblem->grains;
    reinitialiseMatrice(theProblem);
    femCouetteAssemble(theProblem);
    double* B   = theProblem->system->B;
    double* B2  = theProblem->system2->B;

    int i;    
    int n             = myGrains->n;
    double radiusOut  = myGrains->radiusOut;
    double radiusIn   = myGrains->radiusIn;
    
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *norm         = myGrains->norm;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];

// 
// -1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
//
    findElement(theProblem);
    int elem, map[3];
    double xloc[3], yloc[3];
    double vx_loc, vy_loc;
    double iso[2];
    double phi[3];

    for(i = 0; i < n; i++) {
      elem = myGrains->elem[i];
      if(elem == -1){
        vx_loc = 0.0; vy_loc = 0.0;
      }
      else {
        femMeshLocal(theProblem->mesh, elem, map, xloc, yloc);
        isomorphisme(myGrains->x[i], myGrains->y[i], xloc, yloc, iso);
        vx_loc = 0;
        vy_loc = 0;
        double phi1 = 1 - iso[0] - iso[1];
        double phi2 = iso[0];
        double phi3 = iso[1];
        int j;
        vx_loc = phi1 * B[map[0]]  + phi2 * B[map[1]]  + phi3 * B[map[2]];
        vy_loc = phi1 * B2[map[0]] + phi2 * B2[map[1]] + phi3 * B2[map[2]];
      }

      double fx   = m[i] * gx - gamma * vx[i] + gamma * vx_loc;
      double fy   = m[i] * gy - gamma * vy[i] + gamma * vy_loc;
      vx[i]      += fx * dt / m[i];
      vy[i]      += fy * dt / m[i];  
    }

//
// -2- Correction des vitesses pour tenir compte des contacts        
//       

    int iter = 0;
    double error;
           
    do {
        error = femGrainsContactIterate(myGrains,dt,iter);
        iter++; }
    while ((error > tol/dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,error);
    femNorm(vx,vy,norm,n);
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//

    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}

# endif