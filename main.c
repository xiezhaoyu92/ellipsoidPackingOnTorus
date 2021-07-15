//
//  main.c
//  torusPacking
//
//  Created by Zhaoyu Xie on 2/23/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "surface.h"
#include "particles.h"
#include "operation.h"
#include "mt19937ar.h"
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON = 1e-15
#endif

int main(int argc, char * argv[]) {
    
    double R0 = 1; //major radius
    double r0 = 0.5; //minor radius
    //particles
    double a = 0.1; //half length of long axis
    double b = 0.05; //half length of short axis
    int np = 100;
    
    double diffusionTimescale = 0.25;
    
    int animationQ = FALSE;
    int animationRate = 100;
    int readUnfinishedQ = FALSE;
    int readJammedQ = FALSE;
    
    static struct option longOptions[]={
        {"majorAxes",           required_argument,  NULL, 'a'},
        {"minorAxes",           required_argument,  NULL, 'b'},
        {"particleNumber",      required_argument,  NULL, 'n'},
        {"majorRadius",         required_argument,  NULL, 'R'},
        {"minorRadius",         required_argument,  NULL, 'r'},
        {"animate",             no_argument,        NULL, 'm'},
        {"dtDiffusionScale",    required_argument,  NULL, 's'},
        {"readUnfinished",      no_argument,        NULL, 'u'},
        {"readJammed",          no_argument,        NULL, 'j'},
        {0,                     0,                  0,     0 }
    };
    int optIndex = 0;
    int opt;
    while ((opt = getopt_long(argc, argv,"a:b:n:R:r:ms:uj",
                              longOptions, &optIndex )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'a' :
                a = atof(optarg);
                break;
            case 'b' :
                b = atof(optarg);;
                break;
            case 'n' :
                np = atoi(optarg);
                break;
            case 'R' :
                R0 = atof(optarg);
                break;
            case 'r' :
                r0 = atof(optarg);
                break;
            case 'm' :
                animationQ = TRUE;
                break;
            case 's':
                diffusionTimescale = atof(optarg);
                break;
            case 'u' :
                readUnfinishedQ = TRUE;
                break;
            case 'j' :
                readJammedQ = TRUE;
                break;
            default:
                exit(1);
        }
    }
    
    if (readUnfinishedQ && readJammedQ) {
        printf("readJammmed and readUnfinished are mutually exclusive options. Please specify only one.\n");
        exit(1);
    }
    
    double R;
    double r;
    double ROld;
    double rOld;
    
    //undoOverlap
    int jammed = FALSE;
    int maxUndoStepsStart = 1e6;
    double dtOverlap = 1e-4;
    double dtOverlapStartup = 1e-4;
    int maxUndoSteps = 1e4;
    
    double simTimeStart = 0;
    double nextRelaxationTime = 0;
    int relaxationStep = 0;
    double simTime;
    int simStep = 0;
    
    if(readUnfinishedQ||readJammedQ){
        FILE *abFile=NULL;
        abFile = fopen("particleParams.dat","r");
        if (abFile) {
            fscanf(abFile, "%lf %lf", &a, &b);
            fclose(abFile);
        } else {
            printf("abFile pointer is null\n");
            exit(1);
        }
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","r");
        if (npFile) {
            fscanf(npFile, "%i", &np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *rInFile=NULL;
        rInFile = fopen("paramInit.dat", "r");
        if (rInFile) {
            fscanf(rInFile, "%lf %lf", &R0, &r0);
            fclose(rInFile);
        } else {
            printf("rInFile pointer is null\n");
            exit(1);
        }
        if(readUnfinishedQ) {
            FILE *rFileTemp=NULL;
            rFileTemp = fopen("torusParamsTemp.dat","r");
            if (rFileTemp) {
                fscanf(rFileTemp, "%lf %lf", &R, &r);
                fclose(rFileTemp);
            } else {
                printf("rFile pointer is null\n");
                exit(1);
            }
            FILE *simTimeTempFile=NULL;
            simTimeTempFile = fopen("simTimeTemp.dat","r");
            if (simTimeTempFile) {
                fscanf(simTimeTempFile, "%lf", &simTimeStart);
                fclose(simTimeTempFile);
            } else {
                printf("simTimeTempFile pointer is null\n");
                exit(1);
            }
            FILE *nextRelaxationFile=NULL;
            nextRelaxationFile = fopen("nextRelaxationTime.dat","r");
            if (nextRelaxationFile) {
                fscanf(nextRelaxationFile, "%lf", &nextRelaxationTime);
                fclose(nextRelaxationFile);
            } else {
                printf("nextRelaxationFile pointer is null\n");
                exit(1);
            }
            FILE *steps = NULL;
            steps = fopen("stepsTemp.dat","r");
            if(steps) {
                fscanf(steps, "%d %d", &simStep, &relaxationStep);
                fclose(steps);
            }
            else {
                printf("steps pointer is null\n");
                exit(1);
            }
        }
        else {
            FILE *rFile=NULL;
            rFile = fopen("torusParams.dat","r");
            if (rFile) {
                fscanf(rFile, "%lf %lf", &R, &r);
                fclose(rFile);
            } else {
                printf("rFile pointer is null\n");
                exit(1);
            }
            FILE *simTimeArrestedFile=NULL;
            simTimeArrestedFile = fopen("simTimeArrested.dat","r");
            if (simTimeArrestedFile) {
                fscanf(simTimeArrestedFile, "%lf", &simTimeStart);
                fclose(simTimeArrestedFile);
            } else {
                printf("simTimeArrestedFile pointer is null\n");
                exit(1);
            }
            nextRelaxationTime = simTimeStart;
        }
    }
    else {
        FILE *abFile=NULL;
        abFile = fopen("particleParams.dat","w");
        if (abFile) {
            fprintf(abFile, "%.15lf %.15lf\n", a, b);
            fclose(abFile);
        } else {
            printf("abFile pointer is null\n");
            exit(1);
        }
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","w");
        if (npFile) {
            fprintf(npFile, "%i\n", np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *rInFile=NULL;
        rInFile = fopen("paramInit.dat", "w");
        if (rInFile) {
            fprintf(rInFile, "%.15lf %.15lf\n", R0, r0);
            fclose(rInFile);
        } else {
            printf("rInFile pointer is null\n");
            exit(1);
        }
        
        R = R0;
        r = r0;
    }

    particle p[np];
    
    //relaxation parameters
    double Da = pow(2*a/1000, 2)/2; //diffusion coefficient along long axis
    double Db = pow(2*b/1000, 2)/2; //diffusion coefficient along short axis
    double Dbar = (Da+Db)/2; //average diffusion coefficient
    double Dth = pow(2*PI/1000, 2)/2; //diffusion coefficient for angle theta
    //Here 2a 2b and 2pi are the total diffusion without impedement of particles
    double simTimeFinal = 2*a*a/(Da*diffusionTimescale);
    double relaxationSteps = 1e5;
    double dtRelaxation = simTimeFinal/relaxationSteps;
    double dtTol = 1e-5*dtRelaxation;
    int outputSteps = 10000;
    
    double dtDiffusion = 1;
    //double dtDiffusionMin = 1e-3;
    double timestepScaleFactor = 100; // how much to scale the preliminary timestep, which gives an average step equal to the interparticle spacing and results in an acceptance ratio of ~0.5
    double minSpacing = 1e-7; // spacing on which to base minimum dtDiffusion
    double dtDiffusionMin = timestepScaleFactor*minSpacing*minSpacing/(2*Dbar);
    
    if (dtDiffusion>dtRelaxation) {
        printf("Fast relaxation rate: using smaller timestep.\n");
        dtDiffusion = dtRelaxation;
        printf("%e\n",dtDiffusion);
    }
    
    if (readUnfinishedQ){
        FILE *dtFile=NULL;
        dtFile = fopen("dtTemp.dat","r");
        if (dtFile) {
            fscanf(dtFile, "%lf %lf %lf", &dtDiffusion, &dtRelaxation, &dtOverlap);
            fclose(dtFile);
        } else {
            printf("dtFile pointer is null\n");
            exit(1);
        }
    }

    double simTimeLastRelax = simTimeStart;
    double nextRelaxationTimeLastRelax = nextRelaxationTime;
    
    //partition
    int nPart[3];
    nPart[0] = (int)((R0+r0)/a)+1;
    nPart[1] = (int)((R0+r0)/a)+1;
    nPart[2] = (int)(r0/a)+1;
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    double partBoundZ[nPart[2]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = -(R0+r0)+i*2*a;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = -(R0+r0)+i*2*a;
    for(int i=0; i<nPart[2]; i++)
        partBoundZ[i] = -r0+i*2*a;

    clock_t begin, end;
    double timeSpent;
    randinitialize();
    begin = clock();
    
    for (int i=0; i<np; i++) {
        p[i].aa = a;
        p[i].bb = b;
    }

    if(readUnfinishedQ){
        FILE *configurationTemp=NULL;
        configurationTemp = fopen("configurationTemp.asc","r");
        if (configurationTemp) {
            for (int j = 0; j < np; j++) {
                fscanf(configurationTemp, "%lf %lf %lf %lf %lf %lf %lf", &(p[j].position[0]), &(p[j].position[1]), &(p[j].position[2]), &(p[j].director[0]), &(p[j].director[1]), &(p[j].director[2]), &(p[j].position[3]));
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, partBoundZ, nPart);
            }
            fclose(configurationTemp);
        }
        else {
            printf("configuration pointer is null\n");
            exit(1);
        }
    }
    else if(readJammedQ) {
        FILE *configuration=NULL;
        configuration = fopen("configuration.asc","r");
        if (configuration) {
            for (int j = 0; j < np; j++) {
                fscanf(configuration, "%lf %lf %lf %lf %lf %lf %lf", &(p[j].position[0]), &(p[j].position[1]), &(p[j].position[2]), &(p[j].director[0]), &(p[j].director[1]), &(p[j].director[2]), &(p[j].position[3]));
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, partBoundZ, nPart);
            }
            fclose(configuration);
        }
        else {
            printf("configuration pointer is null\n");
            exit(1);
        }
        for(int j=0; j<np; j++) {
            vectorCopy(p[j].position, p[j].postPreRelaxPosition, 4);
            vectorCopy(p[j].director, p[j].postPreRelaxDirector, 3);
            for (int l=0; l<3; l++)
                p[j].postPreRelaxCoord[l] = p[j].coord[l];
        }
    }
    else {
        FILE *init1=NULL; // after projection
        init1 = fopen("init1.asc","w");
        if (!init1) {
            printf("file init1.asc failed to open");
            exit(1);
        }
        
        for(int i=0; i<np; i++) {
            //double theta = 2*PI*genrand_real2();
            //double phi = 2*PI*genrand_real2();
            double theta, phi;
            while(1) {
                phi = 2*PI*genrand_real2();
                double x = genrand_real1();
                if (x<(R+r*cos(phi))/(2*PI*R))
                    break;
            }
            theta = 2*PI*genrand_real2();
            p[i].position[0] = (R+r*cos(phi))*cos(theta);
            p[i].position[1] = (R+r*cos(phi))*sin(theta);
            p[i].position[2] = r*sin(phi);
            p[i].position[3] = PI*genrand_real2();
            
            directorCalc(R, r, p[i].position, p[i].director);
            fprintf(init1, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2], p[i].director[0], p[i].director[1], p[i].director[2], p[i].position[3]);
            
            partitionOneParticle(&(p[i]), partBoundX, partBoundY, partBoundZ, nPart);
        }
        if (init1) fclose(init1);
/*
        FILE *configurationTemp=NULL;
        configurationTemp = fopen("init2S.asc","r");
        if (configurationTemp) {
            for (int j = 0; j < np; j++) {
                fscanf(configurationTemp, "%lf %lf %lf %lf %lf %lf %lf", &(p[j].position[0]), &(p[j].position[1]), &(p[j].position[2]), &(p[j].director[0]), &(p[j].director[1]), &(p[j].director[2]), &(p[j].position[3]));
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, partBoundZ, nPart);
            }
            fclose(configurationTemp);
        }
        else {
            printf("configuration pointer is null\n");
            exit(1);
        }
        printf("result: %d\n", overlapQ(&p[6], &p[158]));
        int i=6;
        int j=158;
        int result=(p[j].coord[0]-p[i].coord[0]==-1||p[j].coord[0]-p[i].coord[0]==0||p[j].coord[0]-p[i].coord[0]==1)&&(p[j].coord[1]-p[i].coord[1]==-1||p[j].coord[1]-p[i].coord[1]==0||p[j].coord[1]-p[i].coord[1]==1)&&(p[j].coord[2]-p[i].coord[2]==-1||p[j].coord[2]-p[i].coord[2]==0||p[j].coord[2]-p[i].coord[2]==1)&&overlapQ(&(p[i]), &(p[j]));
        result = (p[j].coord[2]-p[i].coord[2]==-1||p[j].coord[2]-p[i].coord[2]==0||p[j].coord[2]-p[i].coord[2]==1);
        printf("%d, %d, %d\n",nPart[0],nPart[1],nPart[2]);
        //printf("%d %d\n",p[i].coord[2],p[j].coord[2]);
        printf("result: %d\n",result);
*/        
        undoOverlaps(p, maxUndoStepsStart, dtOverlap, &jammed, a, b, R, r, np, partBoundX, partBoundY,  partBoundZ, nPart);
        if (jammed) {
            printf("overcrowded starting point - use fewer particles or smaller radius\n");
            exit(1);
        }
        
        FILE *init2=NULL; // after relaxation
        init2 = fopen("init2.asc","w");
        if (!init2) {
            printf("file init2.asc failed to open");
            exit(1);
        }
        for (int i = 0; i < np; i++) {
            fprintf(init2, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2], p[i].director[0], p[i].director[1], p[i].director[2], p[i].position[3]);
            
        }
        fclose(init2);
        
        FILE *init3=NULL; // after relaxation
        init3 = fopen("init3.asc","w");
        if (!init3) {
            printf("file init3.asc failed to open");
            exit(1);
        }
        
        //diffuse particles on the surface for startupDiffusionSteps
        int startupDiffusionSteps = 20000;
        double DaStartup = pow(2*a/100, 2)/2;
        double DbStartup = pow(2*b/100, 2)/2;
        double DthStartup = pow(2*PI/100, 2)/2;
        for (int i=0; i<startupDiffusionSteps; i++) {
            if(i%1000==0) printf("startup shuffle: %i\n",i);
            for (int j=0; j<np; j++) {
                int k = randInteger(np)-1;
                vectorCopy(p[k].position, p[k].oldPosition, 4);
                vectorCopy(p[k].director, p[k].oldDirector, 3);
                for (int l=0; l<3; l++)
                    p[k].oldCoord[l] = p[k].coord[l];
                addStochasticForce(&(p[k]), DaStartup, DbStartup, DthStartup, R, r);
                gradDescentConstrainedIntegrate(&(p[k]), dtDiffusion, R, r);
                partitionOneParticle(&(p[k]), partBoundX, partBoundY, partBoundZ, nPart);
                if(anyOverlapQ(k, p, np)||isnan(p[k].position[0])) {
                    vectorCopy(p[k].oldPosition, p[k].position, 4);
                    vectorCopy(p[k].oldDirector, p[k].director, 3);
                    for (int l=0; l<3; l++)
                        p[k].coord[l] = p[k].oldCoord[l];
                }
            }
            /*if (i%animationRate==0)
             for (int j = 0; j < np; j++) {
             fprintf(init3, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2], p[j].director[0], p[j].director[1], p[j].director[2], p[j].position[3]);
             
             }*/
        }
        
        for (int j = 0; j < np; j++) {
            fprintf(init3, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2], p[j].director[0], p[j].director[1], p[j].director[2], p[j].position[3]);
        }
        fclose(init3);
    }

    FILE *animationFile=NULL;
    if (animationQ) {
        animationFile = fopen("animation.dat", "a");
        if (animationFile) {
            fprintf(animationFile, "%.15lf %.15lf\n", R, r);
            for (int j=0; j<np; j++)
                fprintf(animationFile, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2], p[j].director[0], p[j].director[1], p[j].director[2], p[j].position[3]);
            fclose(animationFile);
        }
        else
            printf("animationFile pointer is null\n");
        
    }
    
    FILE* acceptFile = NULL; //record the acceptance rate and dtDiffusion
    int acceptRecordQ = TRUE;
    
    for (simTime = simTimeStart; simTime <= simTimeFinal && !jammed; simTime += dtDiffusion) {
        if (simTime > nextRelaxationTime) {
            ROld = R;
            rOld = r;
            for(int j=0; j<np; j++) {
                vectorCopy(p[j].position, p[j].oldPosition, 4);
                vectorCopy(p[j].director, p[j].oldDirector, 3);
                for (int l=0; l<3; l++)
                    p[j].oldCoord[l] = p[j].coord[l];
            }
            R = R0*(1-simTime/simTimeFinal);
            r = r0*(1-simTime/simTimeFinal);
            for(int j=0; j<np; j++) {
                projectOntoNewSurface(&(p[j]), ROld, rOld, R, r);
                if(isnan(p[j].position[0])) {
                    vectorCopy(p[j].oldPosition, p[j].position, 4);
                    vectorCopy(p[j].oldDirector, p[j].director, 3);
                    for (int l=0; l<3; l++)
                        p[j].coord[l] = p[j].oldCoord[l];
                    projectOntoSurface(&(p[j]), R, r);
                }
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, partBoundZ, nPart);
            }
            int nearJammed = FALSE;
            undoOverlaps(p, maxUndoSteps, dtOverlap, &nearJammed, a, b, R, r, np, partBoundX, partBoundY, partBoundZ, nPart);
            if(nearJammed) {
                R = ROld;
                r = rOld;
                for(int j=0; j<np; j++) {
                    vectorCopy(p[j].postPreRelaxPosition, p[j].position, 4);
                    vectorCopy(p[j].postPreRelaxDirector, p[j].director, 3);
                    for (int l=0; l<3; l++)
                        p[j].coord[l] = p[j].postPreRelaxCoord[l];
                }
                simTime = simTimeLastRelax;
                nextRelaxationTime = nextRelaxationTimeLastRelax;
                dtRelaxation *= 0.5;
                dtOverlap *= 0.5;
                printf("reducing timesteps %e %e\n",dtRelaxation,dtOverlap);
                if (dtRelaxation < dtTol)
                    jammed = TRUE;
                if (dtDiffusion > dtRelaxation)
                    dtDiffusion = dtRelaxation;
            }
            else {
                simTimeLastRelax = simTime;
                nextRelaxationTimeLastRelax = nextRelaxationTime;
                for(int j=0; j<np; j++) {
                    vectorCopy(p[j].position, p[j].postPreRelaxPosition, 4);
                    vectorCopy(p[j].director, p[j].postPreRelaxDirector, 3);
                    for (int l=0; l<3; l++)
                        p[j].postPreRelaxCoord[l] = p[j].coord[l];
                }
                relaxationStep++;
                if (animationQ && relaxationStep % animationRate == 0) {
                    animationFile = fopen("animation.dat", "a");
                    if (animationFile) {
                        fprintf(animationFile, "%.15lf %.15lf\n", R, r);
                        for (int j=0; j<np; j++)
                            fprintf(animationFile, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2], p[j].director[0], p[j].director[1], p[j].director[2], p[j].position[3]);
                        fclose(animationFile);
                    } else
                        printf("animationFile pointer is null\n");
                }
            }
            nextRelaxationTime += dtRelaxation;
            
        }
        //***************************//
        // one set of particle moves //
        //***************************//
        int stepCount = 0;
        int overlapCount = 0;
        for(int i=0; i<np; i++)
            for(int j=0; j<4; j++){
                p[i].force[j] = 0;
                p[i].forceStoc[j] = 0;
            }
        for(int i=0; i<np; i++)
            addStochasticForce(&(p[i]), Da, Db, Dth, R, r);
        
        // going to loop through particles in a random order - set up shuffled array of indices
        int randIndices[np];
        for (int k=0; k<np; k++) {
            randIndices[k]=k;
        }
        for (int k=np-1; k>0; k--) {
            int swapWith = randInteger(k-1);
            int hold = randIndices[k];
            randIndices[k] = randIndices[swapWith];
            randIndices[swapWith] = hold;
        }
        
        for (int k0 = 0; k0 < np; k0++) {
            int k = randIndices[k0];
            vectorCopy(p[k].position, p[k].oldPosition, 4);
            vectorCopy(p[k].director, p[k].oldDirector, 3);
            for (int l=0; l<3; l++)
                p[k].oldCoord[l] = p[k].coord[l];
            gradDescentConstrainedIntegrate(&(p[k]), dtDiffusion, R, r);
            partitionOneParticle(&(p[k]), partBoundX, partBoundY, partBoundZ, nPart);
            stepCount++;
            if(anyOverlapQ(k, p, np)||isnan(p[k].position[0])) {
                vectorCopy(p[k].oldPosition, p[k].position, 4);
                vectorCopy(p[k].oldDirector, p[k].director, 3);
                for (int l=0; l<3; l++)
                    p[k].coord[l] = p[k].oldCoord[l];
                overlapCount++;
            }
        }
        double acceptanceRatio = 1.0 - (double)overlapCount/stepCount;
        
        if (acceptRecordQ) {
            acceptFile = fopen("acceptance.dat","a");
            if (acceptFile) {
                fprintf(acceptFile, "%e %f %e\n", dtDiffusion, acceptanceRatio, dtRelaxation);
                fclose(acceptFile);
            } else {
                printf("acceptFile pointer is null\n");
            }
            
        }
        
        if(acceptanceRatio < 0.5)
            dtDiffusion = dtDiffusion*0.99;
        else
            dtDiffusion = dtDiffusion*1.01;
        if(dtDiffusion < dtDiffusionMin)
            dtDiffusion = dtDiffusionMin;
        if(dtDiffusion > dtRelaxation)
            dtDiffusion = dtRelaxation;
        
        if (simStep%outputSteps==0) {
            FILE *configurationTemp=NULL;
            configurationTemp = fopen("configurationTemp.asc.tmp","w");
            if (configurationTemp) {
                for (int j = 0; j < np; j++)
                    fprintf(configurationTemp, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2], p[j].director[0], p[j].director[1], p[j].director[2], p[j].position[3]);
                fclose(configurationTemp);
            }
            else {
                printf("configuration pointer is null\n");
                exit(1);
            }
            
            FILE *rFileTemp=NULL;
            rFileTemp = fopen("torusParamsTemp.dat.tmp","w");
            if (rFileTemp) {
                fprintf(rFileTemp, "%.15lf %.15lf\n", R, r);
                fclose(rFileTemp);
            } else {
                printf("rFile pointer is null\n");
                exit(1);
            }
            
            FILE *simTimeTempFile=NULL;
            simTimeTempFile = fopen("simTimeTemp.dat.tmp","w");
            if (simTimeTempFile) {
                fprintf(simTimeTempFile, "%lf\n", simTime);
                fclose(simTimeTempFile);
            } else {
                printf("simTimeTempFile pointer is null\n");
                exit(1);
            }
            
            FILE *nextRelaxationFile=NULL;
            nextRelaxationFile = fopen("nextRelaxationTime.dat.tmp","w");
            if (nextRelaxationFile) {
                fprintf(nextRelaxationFile, "%lf\n", nextRelaxationTime);
                fclose(nextRelaxationFile);
            } else {
                printf("nextRelaxationFile pointer is null\n");
                exit(1);
            }
            
            FILE *dtFile=NULL;
            dtFile = fopen("dt.dat.tmp","w");
            if (dtFile) {
                fprintf(dtFile, "%e %e %e\n", dtDiffusion, dtRelaxation, dtOverlap);
                fclose(dtFile);
            } else {
                printf("dtFile pointer is null\n");
                exit(1);
            }
            
            FILE *stepsTemp = NULL;
            stepsTemp = fopen("steps.dat.tmp","w");
            if(stepsTemp) {
                fprintf(stepsTemp, "%d %d\n", simStep+1, relaxationStep);
                fclose(stepsTemp);
            }
            else {
                printf("stepsTemp pointer is null\n");
                exit(1);
            }
            
            
            system("mv -f configurationTemp.asc.tmp configurationTemp.asc");
            system("mv -f torusParamsTemp.dat.tmp torusParamsTemp.dat");
            system("mv -f simTimeTemp.dat.tmp simTimeTemp.dat");
            system("mv -f nextRelaxationTime.dat.tmp nextRelaxationTime.dat");
            system("mv -f dt.dat.tmp dtTemp.dat");
            system("mv -f steps.dat.tmp stepsTemp.dat");
        }
        
        simStep++;
        //printf("simTime/simTimeFinal: %.16lf\n", simTime/simTimeFinal);
    }
    
    end = clock();
    timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    FILE *timeFile=NULL;
    timeFile = fopen("time.dat","w");
    if (timeFile) {
        fprintf(timeFile, "%lf\n", timeSpent);
        fclose(timeFile);
    } else {
        printf("timeFile pointer is null\n");
        exit(1);
    }
    
    FILE *rFile=NULL;
    rFile = fopen("torusParams.dat","w");
    if (rFile) {
        fprintf(rFile, "%.15lf %.15lf\n", R, r);
        fclose(rFile);
    } else {
        printf("rFile pointer is null\n");
        exit(1);
    }
    
    FILE *configuration=NULL;
    configuration = fopen("configuration.asc","w");
    if (configuration) {
        for (int j = 0; j < np; j++)
            fprintf(configuration, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2], p[j].director[0], p[j].director[1], p[j].director[2], p[j].position[3]);
        fclose(configuration);
    }
    else {
        printf("configuration pointer is null\n");
        exit(1);
    }
    
    FILE *simTimeArrestedFile=NULL;
    simTimeArrestedFile = fopen("simTimeArrested.dat","w");
    if (simTimeArrestedFile) {
        fprintf(simTimeArrestedFile, "%lf\n", simTime);
        fclose(simTimeArrestedFile);
    }
    else {
        printf("simTimeArrestedFile pointer is null\n");
        exit(1);
    }
    
    FILE *nextRelaxationArrestedFile=NULL;
    nextRelaxationArrestedFile = fopen("nextRelaxationTimeArrested.dat","w");
    if (nextRelaxationArrestedFile) {
        fprintf(nextRelaxationArrestedFile, "%lf\n", nextRelaxationTime);
        fclose(nextRelaxationArrestedFile);
    }
    else {
        printf("nextRelaxationArrestedFile pointer is null\n");
        exit(1);
    }
    
    FILE *steps = NULL;
    steps = fopen("steps.dat","w");
    if(steps) {
        fprintf(steps, "%d %d\n", simStep, relaxationStep);
        fclose(steps);
    }
    else {
        printf("steps pointer is null\n");
        exit(1);
    }
    
    printf("Hello, World!\n");
    return 0;
}
