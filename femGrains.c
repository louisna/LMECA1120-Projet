/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGrains *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut)
{
    int i,nContact = n*(n-1)/2;
    
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = n;
    theGrains->radiusIn = radiusIn;
    theGrains->radiusOut = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -9.81;
    theGrains->gamma = 0.5;
    
       
    theGrains->x  = malloc(n*sizeof(double));
    theGrains->y  = malloc(n*sizeof(double));
    theGrains->vx = malloc(n*sizeof(double));
    theGrains->vy = malloc(n*sizeof(double));
    theGrains->r  = malloc(n*sizeof(double));
    theGrains->m  = malloc(n*sizeof(double));       
    theGrains->dvBoundary = malloc(n * sizeof(double));
    theGrains->dvContacts = malloc(nContact * sizeof(double));
   
    for(i = 0; i < n; i++) {
        theGrains->r[i] = r;
        theGrains->m[i] = m;
        theGrains->x[i] = (i%5) * r * 2.5 - 5 * r + 1e-8; 
        theGrains->y[i] = (i/5) * r * 2.5 + 2 * r + radiusIn;
        theGrains->vx[i] = 0.0;
        theGrains->vy[i] = 0.0; 
        theGrains->dvBoundary[i] = 0.0; }
 
    for(i = 0; i < nContact; i++)  
        theGrains->dvContacts[i] = 0.0;

  
    return theGrains;
}

femGrains *femGrainsCreateTiny(double radiusIn, double radiusOut)
{
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = 2;
    theGrains->radiusIn   = radiusIn;
    theGrains->radiusOut  = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -10;
    theGrains->gamma = 0.0;
           
    theGrains->x  = malloc(2 * sizeof(double));
    theGrains->y  = malloc(2 * sizeof(double));
    theGrains->vx = malloc(2 * sizeof(double));
    theGrains->vy = malloc(2 * sizeof(double));
    theGrains->r  = malloc(2 * sizeof(double));
    theGrains->m  = malloc(2 * sizeof(double));       
    theGrains->dvBoundary = malloc(2 * sizeof(double));
    theGrains->dvContacts = malloc(sizeof(double));
   
    theGrains->r[0] = 0.1;
    theGrains->r[1] = 0.1;  
    theGrains->m[0] = 1.0;
    theGrains->m[1] = 1.0; 
    theGrains->x[0] = 0.0; 
    theGrains->x[1] = 0.0;
    theGrains->y[0] = -radiusOut + 0.3; 
    theGrains->y[1] = -radiusOut + 0.1; 
    theGrains->vx[0] = 0.0; 
    theGrains->vx[1] = 0.0;            
    theGrains->vy[0] = 0.0; 
    theGrains->vy[1] = 0.0; 
    theGrains->dvBoundary[0] = 0.0;
    theGrains->dvBoundary[1] = 0.0;

    theGrains->dvContacts[0] = 0.0;
     
    return theGrains;
}

void femGrainsFree(femGrains *theGrains)
{
    free(theGrains->x);
    free(theGrains->y);
    free(theGrains->vx);
    free(theGrains->vy);
    free(theGrains->r);
    free(theGrains->m);
    free(theGrains->dvBoundary);
    free(theGrains->dvContacts);
    free(theGrains);
}

double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}
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
            iContact++; }}

    for(i = 0; i < n; i++) {
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
        vy[i] -= dv * ny; }  
    return error;
}

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
