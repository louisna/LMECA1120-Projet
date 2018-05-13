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
    //int    n         = 200;
    int    n         = 20;
    //double radius    = 0.02;
    double radius    = 0.1;
    double mass      = 0.52;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;    
    double dt        = 1e-1;
    double tEnd      = 60.0;
    double tol       = 1e-6;
    double t         = 0;
    double iterMax   = 100;

    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);
   
    femCouetteProblem* theProblem = femCouetteCreate("../data/meshMedium.txt", theGrains);
    double Bfin[theProblem->system->size];
    
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->system->size);

    femCouetteAssemble(theProblem);
    int k;
    //for(k=0;k<theGrains->n;k++){
    //    printf("Bille %d dans element %d\n", k, theGrains->elem[k]);
    //}

    
    GLFWwindow* window = glfemInit("MECA1120 : Projet EF ");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.1;
    char theMessage[256];
    int boolMesh = 1;

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem->mesh,w,h);
        glfemPlotField(theProblem->mesh,theProblem->norm,boolMesh);
        glColor3f(0.0,0.0,0.0);
        
        sprintf(theMessage, "Max : %.4f", femMax(theProblem->system->B,theProblem->system->size));
        glColor3f(1,0,0); glfemDrawMessage(340,460,theMessage);

        for (i=0 ;i < theGrains->n; i++) {     
            glColor3f(0,0,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); 
        }
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn); 
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);  
            femGrainsUpdate(theProblem,dt,tol,iterMax);
            t += dt;
        }
        //while ( glfwGetTime()-currentTime < theVelocityFactor*10) {
            if (glfwGetKey(window,'R') == GLFW_PRESS) 
                theRunningMode = 1; 
            if (glfwGetKey(window,'S') == GLFW_PRESS) 
                theRunningMode = 0; 
            if (glfwGetKey(window,'B') == GLFW_PRESS)
                boolMesh = 1;
            if (glfwGetKey(window,'N') == GLFW_PRESS)
                boolMesh = 0;
        //}

    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1);
           
               
    glfwTerminate(); 
    femCouetteFree(theProblem);
    exit(EXIT_SUCCESS);
    
}

