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
int main(int argc, char** argv) {
    constexpr double length = 16;
    constexpr int nlayers = 4;
    std::vector<lattice::AtomTypes> atoms(nlayers, lattice::AtomTypes::CARBON);
    lattice::PairPotentialProblem uprob(lattice::ljpotent, length, atoms);
    snowgoose::RandomPointGenerator<double> rg(*(uprob.mBox));
    const int n = uprob.mVarTypes.size();
    double x[n];
    double bestx[n];
    constexpr int npoints = 100;
    double v = std::numeric_limits<double>::max();
    for(int i = 0; i < npoints; i ++) {
        rg.getPoint(x);
        double nv = uprob.mObjectives[0]->func(x);
        /*
        std::cout << nv << " = ";
        std::cout << "f(" << snowgoose::VecUtils::vecPrint(n, x) << ")\n";
         */ 
        if(nv < v) {
            v = nv;
            snowgoose::VecUtils::vecCopy(n, x, bestx);
        }
    }
    std::cout << " best x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << " best v = " << v << "\n";
    
    return 0;
}

