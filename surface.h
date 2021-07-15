//
//  surface.h
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#ifndef surface_h
#define surface_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef PI
#define PI 3.141592653589793
#endif


#include <stdio.h>

double thetaTorus(double R, double r, double *position);
double phiTorus(double R, double r, double *position);
void tangent1(double R, double r, double *position, double *t1);
void tangent2(double R, double r, double *position, double *t2);
void normal(double R, double r, double *position, double *n0);

#endif /* surface_h */
