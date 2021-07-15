//
//  surface.c
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#include <math.h>
#include "surface.h"

double thetaTorus(double R, double r, double *position) {
    double x = position[0];
    double y = position[1];
    double z = position[2];

    if(y>=0)
        return acos(x/sqrt(x*x+y*y));
    else
        return 2*PI-acos(x/sqrt(x*x+y*y));
}

double phiTorus(double R, double r, double *position) {
    double x = position[0];
    double y = position[1];
    double z = position[2];
    
    if(z>=0)
        return acos((sqrt(x*x+y*y)-R)/r);
    else
        return 2*PI-acos((sqrt(x*x+y*y)-R)/r);
}

void tangent1(double R, double r, double *position, double *t1) {
    double theta = thetaTorus(R, r, position);
    t1[0] = -sin(theta);
    t1[1] = cos(theta);
    t1[2] = 0;
}

void tangent2(double R, double r, double *position, double *t2) {
    double theta = thetaTorus(R, r, position);
    double phi = phiTorus(R, r, position);
    t2[0] = -cos(theta)*sin(phi);
    t2[1] = -sin(theta)*sin(phi);
    t2[2] = cos(phi);
}

void normal(double R, double r, double *position, double *n0) {
    double theta = thetaTorus(R, r, position);
    n0[0] = (position[0]-R*cos(theta))/r;
    n0[1] = (position[1]-R*sin(theta))/r;
    n0[2] = position[2]/r;
}
