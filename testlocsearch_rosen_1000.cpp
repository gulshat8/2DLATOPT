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
//#include <pairpotentials.hpp>
#include <methods/rosenbrock/rosenbrockmethod.hpp>
#include <methods/advancedcoordescent/advancedcoordescent.hpp>
#include <funccnt.hpp>
#include "tsofproblem.hpp"
#include <fstream>
/*
 * 
 */

void acdSearch(double& v, double* x, const COMPI::MPProblem<double>& prob) {
    LOCSEARCH::AdvancedCoordinateDescent<double> desc(prob);
    desc.getOptions().mSearchType = LOCSEARCH::AdvancedCoordinateDescent<double>::SearchTypes::NO_DESCENT;
//    desc.getOptions().mDoTracing = true;
    
    std::cout << "ADC before v = " << v << std::endl;
    bool rv = desc.search(x, v);
    std::cout << "ADC after v = " << v << std::endl;
}


void rosenSearch(double& v, double* x, const COMPI::MPProblem<double>& prob) {
    LOCSEARCH::RosenbrockMethod<double> desc(prob);
    desc.getOptions().mHInit = std::vector<double>(prob.mBox->mDim, 1.);
    desc.getOptions().mMaxStepsNumber = 10000;
    desc.getOptions().mMinGrad = 1e-3;
    desc.getOptions().mHLB = 1e-4 * desc.getOptions().mMinGrad;
    
    desc.getOptions().mDoOrt = false;
    desc.getOptions().mDoTracing = true;
    
    std::cout << "Rosenbrock before v = " << v << std::endl;
    bool rv = desc.search(x, v);
    std::cout << "Rosenbrock after v = " << v << std::endl;
}

void search(double& v, double* x, const COMPI::MPProblem<double>& prob) {
    //acdSearch(v, x, prob);
    rosenSearch(v, x, prob);
}

#define MULT 1
int main(int argc, char** argv) {
    constexpr double length = 16*MULT;
    constexpr int nlayers = 4;
    std::vector<lattice::AtomTypes> atoms(nlayers, lattice::AtomTypes::CARBON);
//    lattice::PairPotentialProblem uprob(lattice::ljpotent, length, atoms);
    lattice::TsofPotentialProblem uprob(length, atoms);
    constexpr int npoints = 1000;
    snowgoose::RandomPointGenerator<double> rg(*(uprob.mBox), npoints, 1);
    const int n = uprob.mVarTypes.size();
    double x[n];
    double bestx[n];
    double v = std::numeric_limits<double>::max();
    auto fcnt = std::make_shared<COMPI::FuncCnt<double>>(uprob.mObjectives[0]);
    uprob.mObjectives.pop_back();
    uprob.mObjectives.push_back(fcnt);

    int sNUM = 1000;
    double smin = 1.0e9, smax = -1.0e9, smean = 0.0, sstd = 0.0;

	int iter = 1;
	ofstream myfile;
	ofstream myfile2;
	myfile.open ("rosen_1000_local_min.txt");
	myfile2.open("rosen_1000_conf.txt");
	myfile << "Local minimums \n nv,f\n";
	myfile2 << "Configurations \n f()\n";
    while(rg.getPoint(x)){

        double nv = uprob.mObjectives[0]->func(x);
        fcnt->reset();
        search(nv, x, uprob);

/*        for(int k = 0; k < sNUM; k++) {
           for(int j = 4; j < n; j+=3) x[j] = fmod(x[j]+k*x[2], x[j+1]);
           nv = uprob.mObjectives[0]->func(x);
           std::cout << " x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
           std::cout << " nv = " << nv << "\n";
           if(smin > nv) smin = nv;
           if(smax < nv) smax = nv;
           smean += nv/MULT;
           sstd  += (nv/MULT)*(nv/MULT);
        }*/
        smean /= sNUM;
        sstd = sqrt(sstd / sNUM - smean*smean);

//        std::cout << "took " << fcnt->mCounters.mFuncCalls << " function calls\n";
        
        //std::cout << nv << " = ";
        //std::cout << "f(" << snowgoose::VecUtils::vecPrint(n, x) << ")\n";
        
		myfile << iter << "," << nv << "\n";
		myfile2 << iter << "," << snowgoose::VecUtils::vecPrint(n, x) << "\n";
		
        if (nv < v) {
            v = nv;
            snowgoose::VecUtils::vecCopy(n, x, bestx);
        }
		int nleft = npoints - iter;
		std::cout << "Iteration No " + iter + ". Left: " + nleft + "\n";
		++iter;
    }
    /*std::cout << " best x = " << snowgoose::VecUtils::vecPrint(n, bestx) << "\n";
    std::cout << " best v = " << v << "\n";
	myfile << " best x = " << snowgoose::VecUtils::vecPrint(n, bestx) << "\n";
	myfile << " best v = " << v << "\n";*/

    std::cout << "min v  =" << smin/MULT << "\n"; 
    std::cout << "max v  =" << smax/MULT << "\n"; 
    std::cout << "mean v =" << smean << "\n"; 
    std::cout << "std v  =" << sstd << "\n"; 

	myfile.close();
	myfile2.close();
    return 0;
}

