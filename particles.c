//
//  particles.c
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#include "particles.h"
#include <math.h>
#include "operation.h"
#include "surface.h"

void directorCalc(double R, double r, double *position, double *u) {
    double t1[3];
    double t2[3];
    double u1[3];
    double u2[3];
    tangent1(R, r, position, t1);
    tangent2(R, r, position, t2);
    vectorScale(t1, cos(position[3]), u1);
    vectorScale(t2, sin(position[3]), u2);
    vectorAdd(u1, u2, u);
}

//project particle onto current sphere surface, and update the director u. This is used after moving the particles
void projectOntoSurface(particle *p, double R, double r) {
    double x = p->position[0];
    double y = p->position[1];
    double z = p->position[2];
    double origin[3];
    origin[0] = x/sqrt(x*x+y*y)*R;
    origin[1] = y/sqrt(x*x+y*y)*R;
    origin[2] = 0;
    double pos[3];
    pos[0] = x-origin[0];
    pos[1]= y-origin[1];
    pos[2] = z-origin[2];
    double rr = vectorNorm(pos);
    for(int i=0; i<3; i++)
        p->position[i] = pos[i]/rr*r+origin[i];
    directorCalc(R, r, p->position, p->director);
}

//project particle onto new surface. This is used after relaxing the surface.
void projectOntoNewSurface(particle *p, double ROld, double rOld, double R, double r) {
    double theta = thetaTorus(ROld, rOld, p->position);
    double phi = phiTorus(ROld, rOld, p->position);
    p->position[0] = (R+r*cos(phi))*cos(theta);
    p->position[1] = (R+r*cos(phi))*sin(theta);
    p->position[2] = r*sin(phi);
    directorCalc(R, r, p->position, p->director);
}

//compare contact distance and center-to-center distance to judge overlap.
int overlapQ(particle *p1, particle *p2) {
    double *position1 = p1->position;
    double *position2 = p2->position;
    double *u1 = p1->director;
    double *u2 = p2->director;
    double a = p1->aa;
    double b = p1->bb;
    
    double r[3];
    
    vectorSubtract(position2, position1, r);
    double R = vectorNorm(r);
    vectorNormalize(r, r);
    
    double chi = (a*a-b*b)/(a*a+b*b);
    double n1 = vectorDotProduct(r, u1);
    double n2 = vectorDotProduct(r, u2);
    double nu = vectorDotProduct(u1, u2);
    double contactDistance = 2*b/sqrt(1-chi/2*((n1+n2)*(n1+n2)/(1+chi*nu)+(n1-n2)*(n1-n2)/(1-chi*nu)));
    
    if(R < contactDistance)
        return TRUE;
    else
        return FALSE;
}

//for a specific particle i in all particles p, check if it overlaps with its surrounding particles. We define surrounding particles as those in the neighboring partitions
int anyOverlapQ(int i, particle p[], int np) {
    int overlaps = FALSE;
    for(int j=0; j<np; j++) {
        if(j==i)
            continue;
        if((p[j].coord[0]-p[i].coord[0]==-1||p[j].coord[0]-p[i].coord[0]==0||p[j].coord[0]-p[i].coord[0]==1)&&(p[j].coord[1]-p[i].coord[1]==-1||p[j].coord[1]-p[i].coord[1]==0||p[j].coord[1]-p[i].coord[1]==1)&&(p[j].coord[2]-p[i].coord[2]==-1||p[j].coord[2]-p[i].coord[2]==0||p[j].coord[2]-p[i].coord[2]==1)) {
            if(overlapQ(&(p[i]), &(p[j]))) {
                overlaps = TRUE;
                break;
            }
        }
    }
    if(overlaps)
        return TRUE;
    else
        return FALSE;
}
//undoOverlaps
//define the interaction of two overlapped ellipsoids, based on the Gaussian model potential of Berne and Pechukas.
double gaussianModelPotential(double x1, double y1, double z1, double theta1, double x2, double y2, double z2, double theta2, double a, double b, double R, double r) {
    double epsilon0 = 10; //Control the magnitude of potential
    double position1[4] = {x1, y1, z1, theta1};
    double position2[4] = {x2, y2, z2, theta2};
    double u1[3];
    double u2[3];
    directorCalc(R, r, position1, u1);
    directorCalc(R, r, position2, u2);
    
    double rr[3];
    
    vectorSubtract(position2, position1, rr);
    double rrr = vectorNorm(rr);
    vectorNormalize(rr, rr);
    
    double chi = (a*a-b*b)/(a*a+b*b);
    double n1 = vectorDotProduct(rr, u1);
    double n2 = vectorDotProduct(rr, u2);
    double nu = vectorDotProduct(u1, u2);
    double sigma = 2*b*b/(1-chi/2*((n1+n2)*(n1+n2)/(1+chi*nu)+(n1-n2)*(n1-n2)/(1-chi*nu)));
    double epsilon = epsilon0/sqrt(1-chi*chi*nu*nu);
    return epsilon*exp(-rrr*rrr/sigma);
}
//use the gradient of potential for different coordinates to get the force
void addOverlapForce(particle *p1, particle *p2, double a, double b, double R, double r) {
    double x1 = p1->position[0];
    double y1 = p1->position[1];
    double z1 = p1->position[2];
    double theta1 = p1->position[3];
    double x2 = p2->position[0];
    double y2 = p2->position[1];
    double z2 = p2->position[2];
    double theta2 = p2->position[3];
    
    double delta = 0.000001;
    p1->force[0] = p1->force[0]+(gaussianModelPotential(x1-delta,y1,z1,theta1,x2,y2,z2,theta2,a,b,R,r)-gaussianModelPotential(x1+delta,y1,z1,theta1,x2,y2,z2,theta2,a,b,R,r))/2/delta;
    p1->force[1] = p1->force[1]+(gaussianModelPotential(x1,y1-delta,z1,theta1,x2,y2,z2,theta2,a,b,R,r)-gaussianModelPotential(x1,y1+delta,z1,theta1,x2,y2,z2,theta2,a,b,R,r))/2/delta;
    p1->force[2] = p1->force[2]+(gaussianModelPotential(x1,y1,z1-delta,theta1,x2,y2,z2,theta2,a,b,R,r)-gaussianModelPotential(x1,y1,z1+delta,theta1,x2,y2,z2,theta2,a,b,R,r))/2/delta;
    p1->force[3] = p1->force[3]+10*(gaussianModelPotential(x1,y1,z1,theta1-delta,x2,y2,z2,theta2,a,b,R,r)-gaussianModelPotential(x1,y1,z1,theta1+delta,x2,y2,z2,theta2,a,b,R,r))/2/delta;
    p2->force[0] = p2->force[0]+(gaussianModelPotential(x1,y1,z1,theta1,x2-delta,y2,z2,theta2,a,b,R,r)-gaussianModelPotential(x1,y1,z1,theta1,x2+delta,y2,z2,theta2,a,b,R,r))/2/delta;
    p2->force[1] = p2->force[1]+(gaussianModelPotential(x1,y1,z1,theta1,x2,y2-delta,z2,theta2,a,b,R,r)-gaussianModelPotential(x1,y1,z1,theta1,x2,y2+delta,z2,theta2,a,b,R,r))/2/delta;
    p2->force[2] = p2->force[2]+(gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2-delta,theta2,a,b,R,r)-gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2+delta,theta2,a,b,R,r))/2/delta;
    p2->force[3] = p2->force[3]+10*(gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2,theta2-delta,a,b,R,r)-gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2,theta2+delta,a,b,R,r))/2/delta;
}



//move particle, and project onto the surface. Update the information of particle
void gradDescentConstrainedIntegrate(particle *p, double dt, double R, double r) {
    //move particle
    for(int i=0; i<4; i++)
        p->position[i] = p->position[i] + (p->force[i])*dt + (p->forceStoc[i])*sqrt(dt);
    if(p->position[3]>=PI)
        p->position[3] = p->position[3]-PI;
    else if(p->position[3]<0)
        p->position[3] = p->position[3]+PI;
    
    //project onto sphere
    projectOntoSurface(p, R, r);
}

//undooverlaps by checking particle with those in surrounding partitions
void undoOverlaps(particle p[], int maxUndoSteps, double dtOverlap, int *nearJammed, double a, double b, double R, double r, int np, double partBoundX[], double partBoundY[], double partBoundZ[], int nPart[]) {
    for(int i=0; i<np; i++)
        for(int j=0; j<4; j++)
            p[i].forceStoc[j] = 0;
    
    int undoSteps = 0;
    *nearJammed = FALSE;
    int totalOverlapQ;
    do {
        totalOverlapQ = FALSE;
        undoSteps++;
        //printf("%d\n", undoSteps);
        if(undoSteps>maxUndoSteps) {
            *nearJammed = TRUE;
            break;
        }
        for(int i=0; i<np; i++)
            for(int j=0; j<4; j++)
                p[i].force[j] = 0;
        for(int i=0; i<np; i++)
            for(int j=i+1; j<np; j++)
                if((p[j].coord[0]-p[i].coord[0]==-1||p[j].coord[0]-p[i].coord[0]==0||p[j].coord[0]-p[i].coord[0]==1)&&(p[j].coord[1]-p[i].coord[1]==-1||p[j].coord[1]-p[i].coord[1]==0||p[j].coord[1]-p[i].coord[1]==1)&&(p[j].coord[2]-p[i].coord[2]==-1||p[j].coord[2]-p[i].coord[2]==0||p[j].coord[2]-p[i].coord[2]==1)&&overlapQ(&(p[i]), &(p[j]))) {
                    addOverlapForce(&(p[i]), &(p[j]), a, b, R, r);
                    totalOverlapQ = TRUE;
                }
        //project force onto the tangent plane
        for(int i=0; i<np; i++) {
            double n[3];
            normal(R, r, p[i].position, n);
            double projection = vectorDotProduct(p[i].force, n);
            for(int j=0; j<3; j++)
                p[i].force[j] = p[i].force[j] - projection*n[j];
        }

        for(int i=0; i<np; i++) {
            vectorCopy(p[i].position, p[i].oldPosition, 4);
            vectorCopy(p[i].director, p[i].oldDirector, 3);
            gradDescentConstrainedIntegrate(&(p[i]), dtOverlap, R, r);
            if(isnan(p[i].position[0])){
                vectorCopy(p[i].oldPosition, p[i].position, 4);
                vectorCopy(p[i].oldDirector, p[i].director, 3);
            }
        }
        for(int i=0; i<np; i++)
            partitionOneParticle(&(p[i]), partBoundX, partBoundY, partBoundZ, nPart);
        
    } while(totalOverlapQ);
    //printf("undoOverlapSteps are %d\n", undoSteps);
}

//Add stochastic force
void addStochasticForce(particle *p, double Da, double Db, double Dth, double R, double r) {
    double theta = p->position[3];
    double t1[3];
    double t2[3];
    tangent1(R, r, p->position, t1);
    tangent2(R, r, p->position, t2);
    double ETAa, ETAb, ETAth;
    ETAa = randNormal();
    ETAb = randNormal();
    ETAth = randNormal();
    
    for(int i=0; i<3; i++)
        p->forceStoc[i] = (ETAa*cos(theta)*sqrt(2*Da)-ETAb*sin(theta)*sqrt(2*Db))*t1[i] + (ETAa*sin(theta)*sqrt(2*Da)+ETAb*cos(theta)*sqrt(2*Db))*t2[i];
    p->forceStoc[3] = ETAth*sqrt(2*Dth);
}

//Put one particle into its related partition, setting the value of coord[3]
void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], double partBoundZ[], int nPart[]) {
    int i=0, j=0, k=0;
    while(i<nPart[0]&&p->position[0]>=partBoundX[i])
        i++;
    p->coord[0] = i-1;
    
    while(j<nPart[1]&&p->position[1]>=partBoundY[j])
        j++;
    p->coord[1] = j-1;
    
    while(k<nPart[2]&&p->position[2]>=partBoundZ[k])
        k++;
    p->coord[2] = k-1;
}

