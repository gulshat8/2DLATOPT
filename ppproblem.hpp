/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pairpotproblemfact.hpp
 * Author: mposypkin
 *
 * Pair potential problem factory
 * 
 * Created on June 5, 2017, 2:30 PM
 */

#ifndef UNIFORMPP_HPP
#define UNIFORMPP_HPP

#include <vector>
#include <mpproblem.hpp>
#include "atoms.hpp"
#include "latticedata.hpp"
#include "energyfunc.hpp"
#include "pairpotentials.hpp"
#include "latticeutils.hpp"
#include "ppenergy.hpp"

namespace lattice {

    /**
     * The fragment of a lattice of uniform atoms interacting with the pair potential
     */
    class PairPotentialProblem : public COMPI::MPProblem<double> {
    public:

        static void stdBoxFill(const LatticeData& latticeData, snowgoose::Box<double>& box) {
            for (int i = 0; i < latticeData.mNumLayers; i++) {
                const int j = LatticeData::ParamInd::RecSize * i;
                // height
                box.mA[j + LatticeData::ParamInd::Height] = 0.5;
                box.mB[j + LatticeData::ParamInd::Height] = 1;
                // displacement
                if (i == 0) {
                    box.mA[j + LatticeData::ParamInd::Displ] = 0;
                    box.mB[j + LatticeData::ParamInd::Displ] = 0;
                } else {
                    box.mA[j + LatticeData::ParamInd::Displ] = 0;
                    box.mB[j + LatticeData::ParamInd::Displ] = 4;

                }
                // stride
                box.mA[j + LatticeData::ParamInd::Stride] = 0.4;
                box.mB[j + LatticeData::ParamInd::Stride] = 4;
            }
        }

        /**
         * Constructor 
         * @param potent pair potential 
         * @param length the length of the material
         * @param atomTypes atom types (one per layer)
         */
        PairPotentialProblem(PairPotential potent, double length, 
                             const std::vector<lattice::AtomTypes>& atomTypes, 
                             std::function<void(const LatticeData& , snowgoose::Box<double>&) > boxFill = stdBoxFill) {
            const int nlayers = atomTypes.size();
            LatticeData lattice(nlayers);
            lattice.mNumLayers = nlayers;
            lattice.mRadius = lattice::stdCutoff;
            lattice.mLayersAtoms = atomTypes;
            const double qcut = lattice.mRadius * lattice.mRadius;
            const double d = lattice::stdSmothing;
            const double qmin = (lattice.mRadius - d) * (lattice.mRadius - d);
            const double qdelta = qcut - qmin;
            PotentialCutter pc(qcut, qdelta, potent);
            LatticeUtils lut(lattice);
            auto enrg = std::make_shared<PairPotentialEnergy>(lut, pc, length);
            mObjectives.push_back(std::make_shared<EnergyFunc>(enrg));
            const int n = lattice.mNumLayers * LatticeData::ParamInd::RecSize;
            mVarTypes.assign(n, COMPI::MPProblem<double>::VariableTypes::GENERIC);
            mBox = new snowgoose::Box<double>(n);
            boxFill(lattice, *mBox);
        }

    private:

    };

}


#endif /* UNIFORMPP_HPP */

