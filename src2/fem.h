
/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1

typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct 
{
    double *R;
    double *D;
    double *S;
    double *X; 
    double error;      
    int size;
    int iter;        
} femIterativeSolver;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule;
    femIterativeSolver *solverU;
    femIterativeSolver *solverV;
    int size;
    int *number;
    double *soluceU;
    //double *soluceV;
    //double *soluceNorm;
} femDiffusionProblem;


typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule;
    femFullSystem *system;
    femFullSystem *system2;
} femPoissonProblem;

typedef struct {
    int n;
    double radiusIn;
    double radiusOut;
    double gravity[2];
    double gamma;
    double *x;
    double *y;
    double *vx;
    double *vy;
    double *r;
    double *m;
    double *dvBoundary;
    double *dvContacts;
    int *elem; 
} femGrains;


femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femMesh             *femMeshRead(const char *filename);
void                 femMeshWrite(const femMesh* myMesh, const char *filename);
void                 femMeshClean(femMesh *theMesh);
void                 femMeshFree(femMesh *theMesh);
void                 femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y);

femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
void                 femEdgesPrint(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemAlloc(femFullSystem* mySystem, int size);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);
void                 femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femFullSystemGet(femFullSystem* mySystem, int i, int j);

femIterativeSolver*  femIterativeSolverCreate(int size);
void                 femIterativeSolverFree(femIterativeSolver* mySolver);
void                 femIterativeSolverInit(femIterativeSolver* mySolver);
double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
void                 femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue);
void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
int                  femIterativeSolverConverged(femIterativeSolver *mySolver);

femDiffusionProblem *femDiffusionCreate(const char *filename, femRenumType renumType);
void                 femDiffusionFree(femDiffusionProblem *theProblem);
void                 femDiffusionMeshLocal(const femDiffusionProblem *theProblem, const int i, int *map, double *x, double *y, double *u);//, double *v);
void                 femDiffusionCompute(femDiffusionProblem *theProblem);
void                 femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType);

femPoissonProblem   *femPoissonCreate(const char *filename);
void                 femPoissonFree(femPoissonProblem *theProblem);
void                 femPoissonSolve(femPoissonProblem *theProblem);

femGrains  *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut);
femGrains  *femGrainsCreateTiny(double radiusIn, double radiusOut);
void        femGrainsFree(femGrains *myGrains);
void        femGrainsUpdate(femPoissonProblem *theProblem,femGrains *myGrains, double dt, double tol, double iterMax);
double      femGrainsContactIterate(femGrains *myGrains, double dt, int iter);

double               femMin(double *x, int n);
double               femMax(double *x, int n);
void                 femError(char *text, int line, char *file);
void                 femErrorScan(int test, int line, char *file);
void                 femWarning(char *text, int line, char *file);

// Functions of homework 3 :-)



#endif