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

class Cell : public mVector{
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

class Profiles{ // Mathematical Subscripts, index = 1, the first element
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
    double IVAverage(int index, double x1, double x2, int n){ // Initiate values, n = number of points from each interval
        assert(x2 > x1);
        double h = (x2 - x1) / (n - 1);
        double sum = 0;
        switch (index) {
            case 1:
                for (int i = 0; i <= n - 1; i++) {
                    sum += ivDensity(x1 + i * h);
                }
                break;
            case 2:
                for (int i = 0; i <= n - 1; i++) {
                    sum += ivDensity(x1 + i * h) * ivVelocity(x1 + i * h);
                }
                break;
            case 3:
                for (int i = 0; i <= n - 1; i++) {
                    sum += 0.5 * ivDensity(x1 + i * h) * pow(ivVelocity(x1 + i * h), 2)
                    + ivPressure(x1 + i * h) / (gamma - 1);
                }
                break;
            default:
                break;
        }
        return sum / n;
    }
    void InitiateValues(Profiles& values){
        for (int i = 1; i <= nCells; i++) {
            values.u1(i) = IVAverage(1, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
            values.u2(i) = IVAverage(2, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
            values.u3(i) = IVAverage(3, xMin + (i - 1) * xStep, xMin + i * xStep, 10);
        }
    }
    mVector Flux(mVector& u){
        mVector temp(3);
        temp[0] = u[1];
        temp[1] = 0.5 * (3 - gamma) * pow(u[1],2) / u[0] + (gamma - 1) * u[2];
        temp[2] = gamma * u[1] * u[2] / u[0] - 0.5 * (gamma - 1) * pow(u[1],3) / pow(u[0],2);
        return temp;
    }
    mVector GetEigenValues(Cell& cell){
        mVector eigenvalues(3);
        double rho = cell.GetDensity();
        double u = cell.GetVelocity();
        double p = cell.GetPressure();
        double a = sqrt(gamma * p / rho);
        eigenvalues[0] = u;
        eigenvalues[1] = u - a;
        eigenvalues[2] = u + a;
        return eigenvalues;
    }
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
    //------------------LF
    mVector LFFlux(Cell& left, Cell& right, double dt){
        return 0.5 * (Flux(left) + Flux(right)) - (right - left) / (2 * dt / xStep);
    }
    double ComputeTimeStep(Profiles& profile, double atTime){
        double tMin = 0.1;
        for (int i = 1 ; i <= nCells; i++) {
            mVector eigenvalues(3);
            eigenvalues = GetEigenValues(profile[i]);
            double t1 = CFL * xStep / std::abs(eigenvalues[0]);
            double t2 = CFL * xStep / std::abs(eigenvalues[1]);
            double t3 = CFL * xStep / std::abs(eigenvalues[2]);
            if (t1 < tMin) {
                tMin = t1;
            }
            if (t2 < tMin) {
                tMin = t2;
            }
            if (t3 < tMin) {
                tMin = t3;
            }
        }
        if (atTime + tMin > finalTIme) {
            tMin = finalTIme - atTime;
        }
        return tMin;
    }

// ------------------HLL
    mVector HLLFlux(Cell& left, Cell& right, double dt){ // ref Toro
        mVector temp(3);
        mVector LeftEigenvalues = GetEigenValues(left);
        mVector RightEigenvalues = GetEigenValues(right);
        mVector FL = Flux(left);
        mVector FR = Flux(right);
        double SL = min(min(LeftEigenvalues), min(RightEigenvalues));
        double SR = max(max(LeftEigenvalues), max(RightEigenvalues));
        if (SL > 0) {
            temp = FL;
        }
        if (SL < 0 && SR > 0) {
            temp = (SR * FL - SL * FR + SL * SR * (right - left)) / (SR - SL);
        }
        if (SR < 0) {
            temp = FR;
        }
        return temp;
    }

    // ------------------LW
    mVector LWFlux(Cell& left, Cell& right, double dt){
        mVector temp = .5 * (left + right) - 0.5 * dt / xStep * (Flux(right) - Flux(left));
        return Flux(temp);
    }
    // ------------------ FORCE
    mVector FORCEFlux(Cell& left, Cell& right, double dt){
        return 0.5 * (LWFlux(left, right, dt) + LFFlux(left, right, dt));
    }
    void ComputeForward(Profiles& uPre, Profiles& uPost, double dt, mVector (EulerSolver::*pf)(Cell&, Cell&, double)){
        for (int i = 2; i <= nCells - 1; i++) {
            uPost[i] = uPre[i] - dt / xStep * ((this->*pf)(uPre[i], uPre[i + 1], dt) - (this->*pf)(uPre[i - 1], uPre[i], dt));
        }
        uPost[1] = uPre[1] - dt / xStep * ((this->*pf)(uPre[2], uPre[1], dt) - (this->*pf)(uPre[1], uPre[1], dt)); // Profiles has primitive subscripts
        uPost[nCells] = uPre[nCells] - dt / xStep * ((this->*pf)(uPre[nCells], uPre[nCells], dt) - (this->*pf)(uPre[nCells - 1], uPre[nCells], dt));
    }
    void Solve(Profiles& res, mVector (EulerSolver::*method)(Cell&, Cell&, double)){
        ComputeSpatialStep();
        Profiles uPre(nCells);
        Profiles uPost(nCells);
        InitiateValues(uPost);
        double tNow = 0;
        while (startTime + tNow < finalTIme) {
            uPre = uPost;
            double dt = ComputeTimeStep(uPre, tNow);
            ComputeForward(uPre, uPost, dt, method);
            tNow += dt;
        }
        GetOutput(uPost, res);
    }
    void GetOutput(Profiles& uPost, Profiles& res){
        for (int i = 1; i <= nCells; i++) {
            res.u1(i) = GetDensity(uPost, i);
            res.u2(i) = GetVelocity(uPost, i);
            res.u3(i) = GetPressure(uPost, i);
            }
        }
    double min(mVector& vec){
        double min = vec[0];
        for (int i = 1; i != vec.dim(); i++) {
            if(vec[i] < min){
                min = vec[i];
            }
        }
        return min;
    }
    double min(double x1, double x2){
        if (x1 < x2) {
            return x1;
        }
        else {
            return x2;
        }
    }
    double max(double x1, double x2){
        if (x1 > x2) {
            return x1;
        }
        else{
            return x2;
        }
    }
    double max(mVector& vec){
        double max = vec[0];
        for (int i = 1; i != vec.dim(); i++) {
            if (vec[i] > max) {
                max = vec[i];
            }
        }
        return max;
    }
    };


    

#endif /* defined(__EulerEquations__EulerEquations__) */
