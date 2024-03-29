/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"



int main(void)
{   
    int    n = 1;
    double radius    = 0.1;
    double mass      = 0.01;
    double radiusIn  = 0.5;
    double radiusOut = 2.0;    
    double dt      = 1e-1;
    double tEnd    = 8.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100;
    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);
   
    femPoissonProblem* theProblem = femPoissonCreate("../data/meshMedium.txt");
    femRenumType  renumType  = FEM_YNUM;
    femDiffusionProblem* theProblem2 = femDiffusionCreate("../data/meshMedium.txt", renumType);
     
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->system->size);

    
    femDiffusionCompute(theProblem2);
    femPoissonSolve(theProblem);   
    //femSolverPrintInfos(theProblem2->solver);
 
    //printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);
    
    char theMessage[256];
    //sprintf(theMessage, "Max : %.4f", femMax(theProblem->system->B,theProblem->system->size));
    
    GLFWwindow* window = glfemInit("MECA1120 : Projet EF ");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem2->mesh,w,h);
        /*
            faire la norme de soluceU et soluceV
        */
        //printf("U = %f, V= %f; UV = %f \n", theProblem2->soluceU,theProblem2->soluceV,theProblem2->soluceNorm);
        glfemPlotField(theProblem2->mesh,theProblem2->soluceU);
        glColor3f(1.0,0.0,0.0); //glfemDrawMessage(20,460,theMessage);

        for (i=0 ;i < theGrains->n; i++) {     
            glColor3f(1,0,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); }
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn); 
        char theMessage[256];
        //sprintf(theMessage,"Time = %g sec",t);

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            //printf("Time = %4g : ",t);  
            femGrainsUpdate(theProblem,theGrains,dt,tol,iterMax);
            t += dt; }
        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS) 
                theRunningMode = 1; 
          if (glfwGetKey(window,'S') == GLFW_PRESS) 
                theRunningMode = 0; }

    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
           
               
    glfwTerminate(); 
    femPoissonFree(theProblem);
    femDiffusionFree(theProblem2);
    exit(EXIT_SUCCESS);
    
}

