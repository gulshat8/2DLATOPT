/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testuniformpp.cpp
 * Author: mposypkin
 *
 * Created on June 5, 2017, 3:42 PM
 */
#include <iostream>
#include <algorithm>
#include <limits>
#include <pointgen/randpointgen.hpp>
#include <common/vec.hpp>
#include "ppproblem.hpp"
#include "pairpotentials.hpp"

/*
 * 
 */
#define MULT 10000
int main(int argc, char** argv) {
    constexpr double length = 16*MULT;
    constexpr int nlayers = 4;
    std::vector<lattice::AtomTypes> atoms(nlayers, lattice::AtomTypes::CARBON);
    lattice::PairPotentialProblem uprob(lattice::ljpotent, length, atoms);
    snowgoose::RandomPointGenerator<double> rg(*(uprob.mBox), 1000, 201906);
    const int n = uprob.mVarTypes.size();
    double x[n];
    double bestx[n];
    constexpr int npoints = 100;
    double v = std::numeric_limits<double>::max();
    double nv;
    int sNUM = 1000;
    double smin = 1.0e9, smax = -1.0e9, smean = 0.0, sstd = 0.0;
    for(int i = 0; i < 1/*npoints*/; i ++) {
        rg.getPoint(x);
        for(int k = 0; k < sNUM; k++) {
           for(int j = 4; j < n; j+=3) x[j] = fmod(x[j]+k*x[2], x[j+1]);
           nv = uprob.mObjectives[0]->func(x);
           std::cout << " x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
           std::cout << " nv = " << nv << "\n";
           if(smin > nv) smin = nv;
           if(smax < nv) smax = nv;
           smean += nv/MULT;
           sstd  += (nv/MULT)*(nv/MULT);
        }
        smean /= sNUM;
        sstd = sqrt(sstd / sNUM - smean*smean);
        /*
        std::cout << nv << " = ";
        std::cout << "f(" << snowgoose::VecUtils::vecPrint(n, x) << ")\n";
         */ 
        if(nv < v) {
            v = nv;
            snowgoose::VecUtils::vecCopy(n, x, bestx);
        }
    }
//    std::cout << " best x = " << snowgoose::VecUtils::vecPrint(n, bestx) << "\n";
//    std::cout << " best v = " << v << "\n";
    std::cout << "min v  =" << smin/MULT << "\n"; 
    std::cout << "max v  =" << smax/MULT << "\n"; 
    std::cout << "mean v =" << smean << "\n"; 
    std::cout << "std v  =" << sstd << "\n"; 
    return 0;
}

