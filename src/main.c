/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Projet EF
 *  Louis NAVARRE : 1235 16 00
 *  Nahi Michel NASSAR : 1269 16 00
 *
 *  All rights reserved.
 *
 */

#include "glfem.h"



int main(void)
{   
    int    n         = 20;
    double radius    = 0.1;
    double mass      = 0.52;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;    
    double dt        = 0.5e-1;
    double tEnd      = 60.0;
    double tol       = 1e-6;
    double t         = 0;
    double iterMax   = 1000;
    double VEXT      = 4.0;

    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);
   
    femCouetteProblem* theProblem = femCouetteCreate("../data/meshMedium.txt", theGrains);
    theProblem->VEXT = VEXT;
    
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n\n", theProblem->system->size);

    femCouetteAssemble(theProblem);
    
    GLFWwindow* window = glfemInit("MECA1120 : Projet EF ",880,480);
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.00001;
    char theMessage[256];
    int boolMesh = 1, boolPlot = 1;

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem->mesh,w-300,h);

        int height=100,length=780;
        if (boolPlot == 0)
        {
            glfemPlotField(theProblem->mesh,theProblem->system->B,boolMesh);
            sprintf(theMessage,"Max V Grain = %.4f ", femMax(theGrains->vy,n));
            glColor3f(1,0,0); glfemDrawMessage(length,350,theMessage);
            sprintf(theMessage,"Min V Grain = %.4f ", femMin(theGrains->vy,n));
            glColor3f(1,0,0); glfemDrawMessage(length,370,theMessage);
        }
        else if (boolPlot < 0)
        {
            glfemPlotField(theProblem->mesh,theProblem->system2->B,boolMesh);
            sprintf(theMessage,"Max V Grain = %.4f ", femMax(theGrains->vx,n));
            glColor3f(1,0,0); glfemDrawMessage(length,350,theMessage);
            sprintf(theMessage,"Min V Grain = %.4f ", femMin(theGrains->vx,n));
            glColor3f(1,0,0); glfemDrawMessage(length,370,theMessage);
        }
        else 
        {
            glfemPlotField(theProblem->mesh,theProblem->norm,boolMesh);
            sprintf(theMessage,"Max V Grain = %.4f ", femMax(theGrains->norm,n));
            glColor3f(1,0,0); glfemDrawMessage(length,350,theMessage);
            sprintf(theMessage,"Min V Grain = %.4f ", femMin(theGrains->norm,n));
            glColor3f(1,0,0); glfemDrawMessage(length,370,theMessage);
            glColor3f(0.0,0.0,0.0);
        }

        for (i=0 ;i < theGrains->n; i++) {
            double v = theGrains->norm[i];
            glColor3f(v/femMax(theGrains->norm,n),0,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); 
        }

        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);

        sprintf(theMessage, "Max V Fluid : %.4f", femMax(theProblem->system->B,theProblem->system->size));
        glColor3f(1,0,0); glfemDrawMessage(length,400,theMessage);
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);

        sprintf(theMessage,"S - R   = %s - %s ", "Stop","Restart");
        glColor3f(1,0,0); glfemDrawMessage(length,height,theMessage);
        sprintf(theMessage,"O - P   = %s - %s Mesh ", "Active","Desactive");
        glColor3f(1,0,0); glfemDrawMessage(length,height+20,theMessage);
        sprintf(theMessage,"N       = %s Fluid Speed", "Normal");
        glColor3f(1,0,0); glfemDrawMessage(length,height+40,theMessage);
        sprintf(theMessage,"H       = %s Fluid Speed", "Horizontal");
        glColor3f(1,0,0); glfemDrawMessage(length,height+60,theMessage);
        sprintf(theMessage,"V       = %s Fluid Speed", "Vertical");
        glColor3f(1,0,0); glfemDrawMessage(length,height+80,theMessage);
        sprintf(theMessage,"U - D   = %s - %s", "Speed up", "Speed down");
        glColor3f(1,0,0); glfemDrawMessage(length, height+100, theMessage);
        sprintf(theMessage,"Z       = %s", "Set speed to 0.0");
        glColor3f(1,0,0); glfemDrawMessage(length, height+120, theMessage);
        sprintf(theMessage,"Actual speed: %f", theProblem->VEXT);
        glColor3f(1,0,0); glfemDrawMessage(length, height+140, theMessage);


        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);  
            femGrainsUpdate(theProblem,dt,tol,iterMax);
            t += dt;
        }

        if (glfwGetKey(window,'R') == GLFW_PRESS) 
            theRunningMode = 1; 
        if (glfwGetKey(window,'S') == GLFW_PRESS) 
            theRunningMode = 0; 
        if (glfwGetKey(window,'O') == GLFW_PRESS)
            boolMesh = 1;
        if (glfwGetKey(window,'P') == GLFW_PRESS)
            boolMesh = 0;
        if (glfwGetKey(window,'H') == GLFW_PRESS)
            boolPlot = -1;
        if (glfwGetKey(window,'V') == GLFW_PRESS)
            boolPlot = 0;
        if (glfwGetKey(window,'N') == GLFW_PRESS)
            boolPlot = 1;
        if (glfwGetKey(window,'W') == GLFW_PRESS)
            theProblem->VEXT = 0.0;
        if (glfwGetKey(window,'U') == GLFW_PRESS){
            if (theProblem->VEXT >= 10)
                printf("Miximal authorized speed reached\n");
            else
                theProblem->VEXT += 0.4;
        }
        if (glfwGetKey(window,'D') == GLFW_PRESS){
            if (theProblem->VEXT <= -10)
                printf("Minimal authorized speed reached\n");
            else
                theProblem->VEXT -= 0.4;
        }

    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1);
           
               
    glfwTerminate(); 
    femCouetteFree(theProblem);
    exit(EXIT_SUCCESS);
    
}