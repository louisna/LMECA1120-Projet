
#include"fem.h"
/*
Louis NAVARRE : 1235-16-00
Nahi NASSAR : 1269-16-00

Le devoir s'est base en grande partie sur les slides du CM4 (28/2/18).
La structure de femPoissonSolve a été par le professeur au CM5 (7/3/18).
Une ressemblance avec une solution anterieure est possible mais on ne s'est pas base dessus.
*/
#define VEXT 10.0


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
        double value2 = -theMesh->X[theEdges->edges[i].node[j]]*VEXT;;
        femFullSystemConstrain(theSystem, theEdges->edges[i].node[j], value); 
        femFullSystemConstrain(theSystem2, theEdges->edges[i].node[j], value2);     
      }
    }
  }

  //Resolution du systeme par elimination de Gauss
  femFullSystemEliminate(theSystem);
  femFullSystemEliminate(theSystem2);
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
# ifndef NOUPDATE


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
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

    for(i = 0; i < n; i++) {
        double fx = m[i] * gx - gamma * vx[i];
        double fy = m[i] * gy - gamma * vy[i];
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
    printf("iterations = %4d : error = %14.7e \n",iter-1,error);
 
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//

    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}

# endif