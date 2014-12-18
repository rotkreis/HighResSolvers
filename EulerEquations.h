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
#define gamma 1.4

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
public:
    Profiles(int nCells){
        data = std::vector<Cell>(nCells);
        for (int i = 0; i != nCells; i++) {
            data[i] = Cell();
        }
    }
public: // Inheritance
    Cell& operator[](int index){
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
    mVector GetEigenValues(Cell& cell);
    double ComputeTimeStep(Profiles& profile, double atTime);
    
    //------------------ Fluxes --------------------
    mVector Flux(mVector& u);
    mVector LFFlux(Cell& left, Cell& right, double dt);
    mVector HLLFlux(Cell& left, Cell& right, double dt);
    mVector LWFlux(Cell& left, Cell& right, double dt);
    mVector FORCEFlux(Cell& left, Cell& right, double dt);

    void ComputeForward(Profiles& uPre, Profiles& uPost, double dt, mVector (EulerSolver::*method)(Cell&, Cell&, double));
    void Solve(Profiles& res, mVector (EulerSolver::*method)(Cell&, Cell&, double));
    void GetOutput(Profiles& uPost, Profiles& res);
};



    

#endif /* defined(__EulerEquations__EulerEquations__) */
