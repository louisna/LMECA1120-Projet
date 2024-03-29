
#include"fem.h"
/*
Louis NAVARRE : 1235-16-00
Nahi NASSAR : 1269-16-00

Le devoir s'est base en grande partie sur les slides du CM4 (28/2/18).
La structure de femPoissonSolve a été par le professeur au CM5 (7/3/18).
Une ressemblance avec une solution anterieure est possible mais on ne s'est pas base dessus.
*/
#define VEXT 1.0


# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);           
    theProblem->edges = femEdgesCreate(theProblem->mesh);  
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    theProblem->system2 = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femFullSystemFree(theProblem->system2);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
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
    map[j]=theMesh->elem[i*node+j];
    x[j]=theMesh->X[map[j]];
    y[j]=theMesh->Y[map[j]];
  }
}

# endif

# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{

  int option = 1;

  //copie de la structure de theProblem
  femMesh *theMesh=theProblem->mesh;
  femEdges *theEdges=theProblem->edges;
  femDiscrete *theSpace=theProblem->space;
  femIntegration *theRule=theProblem->rule;
  //allocation du systeme (voir slides)
  femFullSystem *theSystem=theProblem->system;
  femFullSystem *theSystem2 = theProblem->system2;
  //terme independant eq Poisson
  double termIndep = 1.0;
  //conditions essentielles homogenes sur la frontiere
  double condHomo = 0.0;

  double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
  int i,j,k,l;
  int map[4];
  /*
  nous construisons un systeme lineaire discret en assemblant
  element par element les matrices A et B de theSystem
  */
  for (i = 0; i < theMesh->nElem; i++)
  {
    //creation d'un maillage local
    femMeshLocal(theMesh,i,map,x,y);

    int nSpace=theSpace->n;

    for (j=0; j < theRule->n; j++)
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
      for (k = 0; k < nSpace; k++)
      {
        /*
        Approximations
        X(xsi,eta)=Sum[x_i*Phi(xsi_i,eta_i)]
        Y(xsi,eta)=Sum[y_i*Phi(xsi_i,eta_i)]
        dX(xsi,eta)=Sum[x_i*dPhi(xsi_i,eta_i)]
        dY(xsi,eta)=Sum[y_i*dPhi(xsi_i,eta_i)]
        (voir slide 15 CM4)
        */
        dxdxsi = dxdxsi + x[k]*dphidxsi[k];  
        dxdeta = dxdeta + x[k]*dphideta[k];   
        dydxsi = dydxsi + y[k]*dphidxsi[k];   
        dydeta = dydeta + y[k]*dphideta[k];
      }

      //calcul du jacobien local (voir slide 17 CM4)
      double J_e = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

      for (k = 0; k < nSpace; k++)
      {
        /*
        derivees fonctions de forme :
            dphidx ne dependent pas de x
            dphidy ne dependent pas de y
        (voir slide 19 CM4)
        */
        dphidx[k] = (dphidxsi[k] * dydeta - dphideta[k] * dydxsi) / J_e;
        dphidy[k] = (dphideta[k] * dxdxsi - dphidxsi[k] * dxdeta) / J_e;
      }

      //remplissage des matrices A et B (voir slides 19-20 CM4) 
      for (k = 0; k < nSpace; k++)
      {
        for(l = 0; l < nSpace; l++)
        {
          double numIntegrateA = (dphidx[k] * dphidx[l] + dphidy[k] * dphidy[l]) * J_e * weight;
          theSystem->A[map[k]][map[l]] = theSystem->A[map[k]][map[l]] + numIntegrateA;
          theSystem2->A[map[k]][map[l]] = theSystem2->A[map[k]][map[l]] + numIntegrateA;;
        }
      }
      for (k = 0; k < nSpace; k++) {
        double numIntegrateB = phi[k] * J_e * weight * termIndep;
        theSystem->B[map[k]] = 0;//theSystem->B[map[k]] + numIntegrateB;
        theSystem2->B[map[k]] = 0;
      }
    }
  }

  //on impose les contraintes sur certains segment
  /*
  for (i=0; i < theEdges->nEdge; i++)
  {
    if (theEdges->edges[i].elem[1] == -1)
    {
      double xe1 = theMesh->X[theEdges->edges[i].node[0]];
      double ye1 = theMesh->Y[theEdges->edges[i].node[0]];
      double xe2 = theMesh->X[theEdges->edges[i].node[1]];
      double ye2 = theMesh->Y[theEdges->edges[i].node[1]];
      double r1 = sqrt(xe1*xe1 + ye1*ye1);
      double r2 = sqrt(xe2*xe2 + ye2*ye2);

      for (j = 0; j < 2; j++)
      {
        //(voir slide 14 CM4)
        if(r1 > 0.5 && r2>0.5){
          femFullSystemConstrain(theSystem, theEdges->edges[i].node[j],VEXT);
        }
        else{
          femFullSystemConstrain(theSystem,theEdges->edges[i].node[j],0);
        }
      }
    }
    */
  int done = 0;
  for(i=0;i<theEdges->nEdge && !done;i++){
    if(theEdges->edges[i].elem[1] == -1){
      for(j=0;j<2;j++){
        double value = theMesh->Y[theEdges->edges[i].node[j]]*VEXT;
        double value2 = -theMesh->X[theEdges->edges[i].node[j]]*VEXT;
        femFullSystemConstrain(theSystem, theEdges->edges[i].node[j], value); 
        femFullSystemConstrain(theSystem2, theEdges->edges[i].node[j], value2);
      }
    }
  }

  //Resolution du systeme par elimination de Gauss
  femFullSystemEliminate(theSystem);
  femFullSystemEliminate(theSystem2);

  //for(i=0;i<theSystem->size;i++){
  //  theSystem->B[i] = sqrt(theSystem->B[i]*theSystem->B[i] + theSystem2->B[i]*theSystem2->B[i]);
  //}
}

# endif

void isomorphisme(double x, double y, double *x_loc, double *y_loc, double *to_return){
  double x_iso, y_iso;
  x_iso = -(x*y_loc[0] - x_loc[0]*y - x*y_loc[2] + x_loc[2]*y + x_loc[0]*y_loc[2] - x_loc[2]*y_loc[0])/(x_loc[0]*y_loc[1] - x_loc[1]*y_loc[0] - x_loc[0]*y_loc[2] + x_loc[2]*y_loc[0] + x_loc[1]*y_loc[2] - x_loc[2]*y_loc[1]);
  y_iso =  (x*y_loc[0] - x_loc[0]*y - x*y_loc[1] + x_loc[1]*y + x_loc[0]*y_loc[1] - x_loc[1]*y_loc[0])/(x_loc[0]*y_loc[1] - x_loc[1]*y_loc[0] - x_loc[0]*y_loc[2] + x_loc[2]*y_loc[0] + x_loc[1]*y_loc[2] - x_loc[2]*y_loc[1]);

  to_return[0] = x_iso;
  to_return[1] = y_iso;
}

# ifndef FINDELEMENT

void findElement(femGrains *theGrains, femPoissonProblem *theProblem){
  femMesh *theMesh = theProblem->mesh;
  int n_g     = theGrains->n;
  int *elem_g = theGrains->elem;
  double *x_g    = theGrains->x;
  double *y_g    = theGrains->y;

  double *X_m    = theMesh->X;
  double *Y_m    = theMesh->Y;
  int *elem_m = theMesh->elem;
  int n_e     = theMesh->nElem;

  int i,j;
  double x_elem[3];
  double y_elem[3];
  int map[4];
  double iso[2];

  for(i=0;i<n_g;i++){
    double x_loc = x_g[i];
    double y_loc = y_g[i];
    int dedans = 0;

    for(j=0;j<n_e && dedans == 0;j++){
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

# ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  
{
    int i,j,iContact;    
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
    
    double error = 0.0;
    double gap,rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;

    error = 0;
    iContact = 0;
    for(i = 0; i < n; i++) {
        for(j = i+1; j < n; j++) { 
            rx = (x[j]-x[i]);
            ry = (y[j]-y[i]);
            rr = sqrt(rx*rx+ry*ry);
            nx = rx/rr;
            ny = ry/rr;
            if (iter == 0) {
                  dv = dvContacts[iContact]; }
            else {
                  vn = (vx[i]-vx[j])*nx + (vy[i]-vy[j])*ny ;
                  gap = rr - (r[i]+r[j]);
                  dv = fmax(0.0, vn + dvContacts[iContact] - gap/dt);
                  dv = dv - dvContacts[iContact];                      
                  dvContacts[iContact] += dv; 
                  error = fmax(fabs(dv),error); }
            vx[i] -= dv * nx * m[j] / ( m[i] + m[j] );
            vy[i] -= dv * ny * m[j] / ( m[i] + m[j] );
            vx[j] += dv * nx * m[i] / ( m[i] + m[j] );
            vy[j] += dv * ny * m[i] / ( m[i] + m[j] );                  
            iContact++; 
        }
        rr = sqrt(x[i]*x[i]+y[i]*y[i]);      
        nx = x[i]/rr;
        ny = y[i]/rr;
        if (iter == 0) {
            dv = dvBoundary[i]; }
        else {
            vn = vx[i]*nx + vy[i]*ny ;      
            gap = rOut - rr - r[i];
            dvOut = fmax(0.0, vn + dvBoundary[i] - gap/dt);
            gap = rr - rIn - r[i];
            dvIn  = fmax(0.0,-vn - dvBoundary[i] - gap/dt);
            dv = dvOut - dvIn - dvBoundary[i]; 
            dvBoundary[i] += dv; 
            error = fmax(fabs(dv),error); }
        vx[i] -= dv * nx;
        vy[i] -= dv * ny;
    }
    return error;
}

# endif

double *theGlobalCoord;

int compare(const void *nodeOne, const void *nodeTwo) 
{
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobalCoord[*iOne] - theGlobalCoord[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
    return  0;  
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i, *inverse;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*theProblem->mesh->nNode);
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theProblem->mesh->X;
            qsort(inverse, theProblem->mesh->nNode, sizeof(int), compare);
            for (i = 0; i < theProblem->mesh->nNode; i++)
                theProblem->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*theProblem->mesh->nNode);
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theProblem->mesh->Y;
            qsort(inverse, theProblem->mesh->nNode, sizeof(int), compare);
            for (i = 0; i < theProblem->mesh->nNode; i++)
                theProblem->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}


void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, 
                                double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j,myRow;
    if (mySolver->iter == 0) {
        for (i = 0; i < nLoc; i++) { 
            myRow = map[i];
            mySolver->R[myRow] -= Bloc[i];
            for(j = 0; j < nLoc; j++) {
                mySolver->R[myRow] += Aloc[i*nLoc+j]*Uloc[j];}}}
    
    for (i = 0; i < nLoc; i++) {
        myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            mySolver->S[myRow] += Aloc[i*nLoc+j] * mySolver->D[myCol]; }}
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) 
{

//
// Pour conserver les conditions essentielles homogenes,
// il suffit d'annuler les deux increments pour le residu et la direction :-)

// Observer qu'il est impossible de resoudre le probleme avec des conditions ess. non-homogenes
// tel que c'etait implemente dans le devoir :-)
// Il faudrait modifier tres legerement la fonction femDiffusionCompute 
// pour que cela soit possible :-)
//

        mySolver->R[myNode] = 0.0;
        mySolver->S[myNode] = 0.0; 
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{


// R   = vecteur R (residu)
// X   = dX
// D   = vecteur D (direction)
// S   = vecteur A D (direction conjugee)

// Demarrage de l'algorithme
// La direction initiale est le residu et on initialise l'increment a zero.

    mySolver->iter++;
    double error = 0.0; int i;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]); }
    
  
    if (mySolver->iter == 1) {
        for (i=0; i < mySolver->size; i++) {
            mySolver->X[i] = 0; 
            mySolver->D[i] = mySolver->R[i]; }}           
    else {    
        double denAlpha = 0.0;
        for (i=0; i < mySolver->size; i++) {
            denAlpha += mySolver->R[i] * mySolver->S[i]; }

        double alpha = -error/denAlpha;
    
        double numBeta = 0.0;
        for (i=0; i < mySolver->size; i++) {            
            mySolver->R[i] = mySolver->R[i] + alpha * mySolver->S[i];
            numBeta += mySolver->R[i] * mySolver->R[i]; }
            
        double beta = numBeta/error;
        for (i=0; i < mySolver->size; i++) {
            mySolver->X[i] = alpha * mySolver->D[i];
            mySolver->D[i] = mySolver->R[i] + beta * mySolver->D[i]; 
            mySolver->S[i] = 0.0; }}
     
    mySolver->error = sqrt(error);
    return(mySolver->X);
}

# ifndef NOUPDATE

void femGrainsUpdate(femPoissonProblem *theProblem, femGrains *myGrains, double dt, double tol, double iterMax)
{

    //femPoissonSolve(theProblem);
    double* B = theProblem->system->B;
    double* B2 = theProblem->system2->B;

    int i;    
    int n = myGrains->n;
    
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];

// 
// -1- Calcul des nouvelles vitesses des grains sur base de la gravit� et de la trainee
//
    findElement(myGrains, theProblem);

    for(i = 0; i < n; i++) {
      //printf("Vitesses : %f %f\n", B[myGrains->elem[i]], B2[myGrains->elem[i]]);
        double fx = m[i] * 0 - gamma * vx[i] + gamma * B[myGrains->elem[i]];
        double fy = m[i] * gy - gamma * vy[i] + gamma * B2[myGrains->elem[i]];
        vx[i] += fx * dt / m[i];
        vy[i] += fy * dt / m[i];  }

//
// -2- Correction des vitesses pour tenir compte des contacts        
//       

    int iter = 0;
    double error;
           
    do {
        error = femGrainsContactIterate(myGrains,dt,iter);
        iter++; }
    while ((error > tol/dt && iter < iterMax) || iter == 1);
    //printf("iterations = %4d : error = %14.7e \n",iter-1,error);
 
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//

    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}

# endif