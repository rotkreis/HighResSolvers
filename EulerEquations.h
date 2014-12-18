//
//  EulerEquations.h
//  EulerEquations
//
//  Created by Li Xinrui on 11/23/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __EulerEquations__EulerEquations__
#define __EulerEquations__EulerEquations__

#include <stdio.h>
#include <cassert>
#include "Matrix.h"

typedef double (*pFunc)(double x);
typedef mVector (*pLimiter)(const mVector& r);

#define gamma 1.4
double min(const mVector& vec);
double max(const mVector& vec);
double min(double x1, double x2);
double max(double x1, double x2);
double absMax(const mVector& v);
double GetDensity(const mVector& vec);
double GetVelocity(const mVector& vec);
double GetPressure(const mVector& vec);

class Cell : public mVector{ // Each Cell
public:
    Cell():
    mVector(3) {}
public: // Inheritance
    using mVector::operator=;
    using mVector::operator*=;
public:
    double GetDensity(){
        return (*this)[0];
    }
    double GetVelocity(){
        return (*this)[1]/(*this)[0];
    }
    double GetPressure(){
        return (gamma - 1) * ((*this)[2] - 0.5 * pow((*this)[1], 2) / (*this)[0]);
    }
};


class Profiles{ // Solution, Mathematical Subscripts, index = 1, the first element
public:
    std::vector<Cell> data;
    int nCells;
public:
    Profiles(int nCells):
    nCells(nCells){
        data = std::vector<Cell>(nCells);
        for (int i = 0; i != nCells; i++) {
            data[i] = Cell();
        }
    }
public: // Inheritance
    // Boundary
    Cell& operator[](int index){
        if (index < 1) {
            return data[0];
        }
        if (index > nCells) {
            return data[nCells - 1];
        }
        return data[index - 1];
    }
    const Cell& operator[](int index) const{
        if (index < 1) {
            return data[0];
        }
        if (index > nCells) {
            return data[nCells -1];
        }
        return data[index - 1];
    }
public:
    double& u1(int index){
        return data[index - 1][0];
    }
    double& u2(int index){
        return data[index - 1][1];
    }
    double& u3(int index){
        return data[index - 1][2];
    }
};

class EulerSolver{
public:
    int nCells;
    double xStep;
    double timeStep;
    pFunc ivDensity;
    pFunc ivPressure;
    pFunc ivVelocity;
    double startTime;
    double finalTIme;
    double xMin;
    double xMax;
    double CFL = 0.8;
    // Constructor
public:
    EulerSolver(pFunc rho, pFunc p, pFunc u):
    ivDensity(rho), ivPressure(p),ivVelocity(u){}
    // Initiate
public:
    void SetRange(double x1, double x2){
        xMin = x1;
        xMax = x2;
    }
    void SetTime(double t1, double t2){
        startTime = t1;
        finalTIme = t2;
    }
    void SetCellNumber(int n){
        nCells = n;
    }
    void SetCFL(double x){
        CFL = x;
    }
    // Preparations
public:
    void ComputeSpatialStep(){
        xStep = (xMax - xMin) / nCells;
    }
    double IVAverage(int index, double x1, double x2, int n);
    void InitiateValues(Profiles& values);
    double GetDensity(Profiles& profile, int index){
        return profile[index].GetDensity();
    }
    double GetVelocity(Profiles& profile, int index) {
        return profile[index].GetVelocity();
    }
    double GetPressure(Profiles& profile, int index){
        return profile[index].GetPressure();
    }
    // Computations
public:
    mVector GetEigenValues(const mVector& cell);
    double ComputeTimeStep(Profiles& profile, double atTime);
    
    //------------------ Fluxes --------------------
    mVector Flux(const mVector& u);
    mVector LFFlux(Cell& left,Cell& right, double dt);
    mVector HLLFlux(Cell& left, Cell& right, double dt);
    mVector LWFlux(Cell& left, Cell& right, double dt);
    mVector FORCEFlux(Cell& left, Cell& right, double dt);

    void ComputeForward(Profiles& uPre, Profiles& uPost, double dt, mVector (EulerSolver::*method)(Cell&, Cell&, double));
    void Solve(Profiles& res, mVector (EulerSolver::*method)(Cell&, Cell&, double));
    void GetOutput(Profiles& uPost, Profiles& res);
    
    // ------------- HighRes --------
    mVector componentWiseMultiply(const mVector& v1, const mVector& v2){
        assert(v1.dim() == v2.dim());
        mVector temp(v1.dim());
        for (int i = 0; i != temp.dim(); i++) {
            temp[i] = v1[1] * v2[i];
        }
        return temp;
    }
    mVector vectorR(const Profiles& u, int index){
        mVector r(3);
        for (int i = 0; i != 3; i++) {
            if (u[index][i] - u[index - 1][i] == 0) {
                r[i] = 0;
            }
            else {
            r[i] = (u[index] - u[index - 1])[i] / (u[index + 1] - u[index])[i];
            }
        }
        return r;
    }
    // hard to name for L R and i +- 1/2! !!!! !!! oh me !!
    mVector uNewLeft(const Profiles& u, int index, pLimiter limiter){// index for i + 1/2, index- 1 for i - 1/2
//        for (int i = 0; i != 3; i++) {
//            assert((u[index+1] - u[index])[i] != 0);
//        }
        mVector r = vectorR(u, index);
        return u[index] + 0.5 * componentWiseMultiply(limiter(r), u[index+1] - u[index]);
    }
    mVector uNewRight(const Profiles& u, int index, pLimiter limiter){ // index for i + 1/2
//        for (int i = 0; i != 3; i++) {
//            assert((u[index+2] - u[index+1])[i] != 0);
//        }
        mVector r = vectorR(u, index + 1);
        return u[index + 1] - 0.5 * componentWiseMultiply(limiter(r), u[index + 2] - u[index - 1]);
    }
    mVector KTNumericalFlux(const Profiles& u, int index, pLimiter limiter){// index for i + 1/2, index - 1 for i - 1/2
        mVector temp(3);
        double a = max(absMax(GetEigenValues(u[index])), absMax(GetEigenValues(u[index + 1])));
        temp = (Flux(uNewRight(u, index, limiter)) + Flux(uNewLeft(u, index, limiter)));
        temp -= a * (uNewRight(u, index, limiter) - uNewLeft(u, index, limiter));
        temp *= 0.5;
        return temp;
    }
    void HighResComputeForward(const Profiles& uPre, Profiles& uPost, double dt, pLimiter limiter){
        for (int i = 1; i <= nCells; i++) {
            uPost[i] = uPre[i] - dt / xStep * (KTNumericalFlux(uPre, i, limiter) - KTNumericalFlux(uPre, i - 1, limiter));
        }
    }
    Profiles HighResSolve(pLimiter limiter){
        ComputeSpatialStep();
        Profiles uPre(nCells);
        Profiles uPost(nCells);
        InitiateValues(uPost);
        double tNow = startTime;
        while (tNow < finalTIme) {
            uPre = uPost;
            double dt = ComputeTimeStep(uPre, tNow);
            HighResComputeForward(uPre, uPost, dt, limiter);
            tNow += dt;
        }
        return uPost;
    }
};

namespace Limiter {
    double minmod(double r);
    mVector minmod(const mVector& R);
    double superbee(double r);
    mVector superbee(const mVector& R);
    double vanLeer(double r);
    mVector vanLeer(const mVector& R);
}


    

#endif /* defined(__EulerEquations__EulerEquations__) */
