//
//  particles.h
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#ifndef particles_h
#define particles_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#include <stdio.h>

typedef struct {
    double position[4]; //x, y, z are the first three terms of position, and angle theta is the fourth.
    double oldPosition[4];
    double director[3];
    double oldDirector[3];
    double postPreRelaxPosition[4];
    double postPreRelaxDirector[3];
    
    double force[4];
    double forceStoc[4];
    double aa;
    double bb;
    int oldCoord[3];
    int coord[3];
    double postPreRelaxCoord[3];
} particle;

void directorCalc(double R, double r, double *position, double *u);
void projectOntoSurface(particle *p, double R, double r);
void projectOntoNewSurface(particle *p, double ROld, double rOld, double R, double r);

int overlapQ(particle *p1, particle *p2);
int anyOverlapQ(int i, particle p[], int np);

double gaussianModelPotential(double x1, double y1, double z1, double theta1, double x2, double y2, double z2, double theta2, double a, double b, double R, double r);
void addOverlapForce(particle *p1, particle *p2, double a, double b, double R, double r);
void gradDescentConstrainedIntegrate(particle *p, double dt, double R, double r);
void undoOverlaps(particle p[], int maxUndoSteps, double dtOverlap, int *nearJammed, double a, double b, double R, double r, int np, double partBoundX[], double partBoundY[], double partBoundZ[], int nPart[]);
void addStochasticForce(particle *p, double Da, double Db, double Dth, double R, double r);

void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], double partBoundZ[], int nPart[]);

#endif /* particles_h */
