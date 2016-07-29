/**************************************************************
 *                                                            *
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' Background Suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2016  Julian Pychy                          *
 *                                                            *
 *                                                            *
 *  Description:                                              *
 *                                                            *
 *  License:                                                  *
 *                                                            *
 *  This file is part of WiBaS                                *
 *                                                            *
 *  WiBaS is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public       *
 *  License as published by the Free Software Foundation,     *
 *  either version 3 of the License, or (at your option) any  *
 *  later version.                                            *
 *                                                            *
 *  WiBaS is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied        *
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR   *
 *  PURPOSE. See the GNU General Public License for more      *
 *  details.                                                  *
 *                                                            *
 *  You should have received a copy of the GNU General        *
 *  Public License along with WiBaS (license.txt). If not,    *
 *  see <http://www.gnu.org/licenses/>.                       *
 *                                                            *
 *************************************************************/


#include <iostream>
#include <cmath>
#include <algorithm>

#include "WibasCore.hh"
#include "FitResult.hh"
#include "PhasespacePoint.hh"
#include "WibFitFunction.hh"
#include "FastPointMap.hh"

#include "RooMsgService.h"



WiBaS::WiBaS(WibFitFunction& pfitFunction) :
    PhasespacePointCloud(1),
    numNearestNeighbors(200),
    fitFunction(&pfitFunction)
{
    RooMsgService::instance().setSilentMode(true);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
}



void WiBaS::SetNearestNeighbors(unsigned int pnumNearestNeighbors){

    numNearestNeighbors = pnumNearestNeighbors;
}



void WiBaS::AddPhasespacePoint(PhasespacePoint& newPhasespacePoint){

    if(!CheckMassInRange(newPhasespacePoint)){
        *_qout << "WARNING: Attempt to add a particle outside mass range (m1="
               << newPhasespacePoint.GetMass() << ", m2="
               << newPhasespacePoint.GetMass2() << "). "
               << "Rejected." << std::endl;
        return;
    }

    PhasespacePointCloud::AddPhasespacePoint(newPhasespacePoint, 1);
}



bool WiBaS::CalcWeight(PhasespacePoint &refPhasespacePoint){

    ArrangePointCoordinates(refPhasespacePoint);

    if(!CheckMassInRange(refPhasespacePoint)){
        *_qout << "WARNING: Attempt to calculate weight of a particle outside mass range (m1="
               << refPhasespacePoint.GetMass() << ", m2="
               << refPhasespacePoint.GetMass2() << "). " << std::endl;
        return false;
    }


    // calculate phasespace distances
    std::vector<FastPointMap> pointMapVector;
    std::vector<PhasespacePoint*>::iterator it;
    std::vector<PhasespacePoint*> phasespacePointVector = GetPointVector();


    for (it=phasespacePointVector.begin(); it!=phasespacePointVector.end(); ++it){
        FastPointMap newMapEntry((*it), CalcPhasespaceDistance((*it), &refPhasespacePoint));
        pointMapVector.push_back(newMapEntry);
    }


    // sort list
    FastPointMap compHelper(NULL, 0);
    std::sort(pointMapVector.begin(), pointMapVector.end(), compHelper);


    // Cut vector
    int cutIndex=-1;
    double weightsum=0;

    for(unsigned int i=1; i<pointMapVector.size();i++){

        weightsum += pointMapVector.at(i)._phasespacePoint->GetInitialWeight();

        if(weightsum >= numNearestNeighbors){
            cutIndex = i;
            break;
        }
    }

    if(cutIndex <= 0){
        *_qout << "ERROR: Too few events available for numNearestNeighbors = " << numNearestNeighbors << std::endl;
        return false;
    }

   pointMapVector.resize(cutIndex + 1);


    // Fill the fit function with the neighbor data
    std::vector<FastPointMap>::iterator it2;
    for(it2=pointMapVector.begin() + 1; it2!=pointMapVector.end(); ++it2) // skip nearest event (=ref event?, TODO: check this!)
    {
        fitFunction->AddData(*((*it2)._phasespacePoint));
    }


    // Do the fit
    FitResult* fitResult = fitFunction->DoFit(refPhasespacePoint.GetMass(),  refPhasespacePoint.GetMass2());


    // Check fit result
    if((fitResult == NULL) || (fitResult->rooFitResult == NULL) || (fitResult->rooFitResult->status() != 0)){
        *_qout << "ERROR: Fit did not converge or returned NULL pointer" << std::endl;
        delete fitResult;
        return false;
    }

    int covQual = fitResult->rooFitResult->covQual();

    if(covQual == 2){
        *_qout << "INFO: covariance matrix forced positive-definite" << std::endl;
    }
    else if(covQual == 1){
        *_qout << "WARNING: covariance matrix not accurate" << std::endl;
    }
    else if(covQual != 3){
        *_qout << "WARNING: covQual = " << covQual << std::endl;
    }


    //Get weight at m0
    double Q = fitResult->weight;

    if(Q > 1) {
        *_qout << "WARNING: Q > 1. Setting Q = 1." << std::endl;
        Q = 1.0;
    }
    else if(Q < 0){
        *_qout << "WARNING: Q < 0. Setting Q = 0." << std::endl;
        Q = 0.0;
    }

    refPhasespacePoint.SetWeight(Q);
    refPhasespacePoint.SetWeightError(fitResult->weightError);

    delete fitResult;
    return true;
}



void WiBaS::SaveNextFitToFile(std::string fileName){

    fitFunction->SaveNextFitToFile(fileName);
}



bool WiBaS::CheckMassInRange(PhasespacePoint &refPhasespacePoint) const {

    double minMass = fitFunction->GetMinMass();
    double maxMass = fitFunction->GetMaxMass();

    if(refPhasespacePoint.GetMass() < minMass || refPhasespacePoint.GetMass() > maxMass){
     return false;
    }

    if(refPhasespacePoint.IsMass2Set() &&
      (refPhasespacePoint.GetMass2() < minMass || refPhasespacePoint.GetMass2() > maxMass)){
        return false;
    }

    return true;
}



void WiBaS::SetCalcErrors(bool set){

    fitFunction->SetCalcErrors(set);
}


