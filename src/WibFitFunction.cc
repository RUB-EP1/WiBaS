/**************************************************************
 *                                                            *
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' Background Suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2015  Julian Pychy                          *
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

#include "WibFitFunction.hh"
#include "FitResult.hh"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"



WibFitFunction::WibFitFunction(double pminMass, double pmaxMass) :
    data(NULL),
    totalIntensity(NULL),
    calcError(false),
    saveNextFitToFile(false),
    minMass(pminMass),
    maxMass(pmaxMass)
{
    mass = new RooRealVar("mass", "mass", 0, maxMass - minMass);
    mass2 = new RooRealVar("mass2", "mass2", 0, maxMass - minMass);
    initialWeight = new RooRealVar("initialWeight", "initialWeight", 0);
}



WibFitFunction::~WibFitFunction(){

    delete mass;
    delete mass2;
    delete initialWeight;

    if(data != NULL)
        delete data;

    if(totalIntensity != NULL)
        delete totalIntensity;

}



bool WibFitFunction::GetCalcError() const {

    return calcError;
}



void WibFitFunction::SetCalcErrors(bool set){

    calcError = set;
}



FitResult* WibFitFunction::DoFit(double eventMass, double eventMass2){

    if(data == NULL)
        return NULL;


    FitResult* fitResult = DoFitD(eventMass - minMass, eventMass2 - minMass);

    if(saveNextFitToFile){
        SaveFitToFile(saveNextFitFileName);
        saveNextFitToFile = false;
    }

    data->reset();

    // Check parameters
    RooFitResult* rooFitResult = fitResult->rooFitResult;
    const RooArgList finalParams = rooFitResult->floatParsFinal();
    const int nFreeParams = finalParams.getSize();
    const RooArgList initialParams = rooFitResult->floatParsInit();
    RooArgList unorderedLocalParams = GetParamList();

    for(int i=0; i<nFreeParams; i++){
        RooRealVar* currentRefVar = dynamic_cast<RooRealVar*>(initialParams.at(i));
        int index = unorderedLocalParams.index(currentRefVar->GetName());

        if(index < 0){
            std::cout << "ERROR: could not find parameter " << currentRefVar->GetName()
                      << " in ::GetParamList() List" << std::endl;
            return NULL;
        }
    }

    // Calculate error
    if(GetCalcError()){
        std::vector<double> derivatives;

        for(int i=0; i<nFreeParams; i++){
            RooRealVar* currentRefVar = dynamic_cast<RooRealVar*>(finalParams.at(i));              // Insert better idea to match
            int index = unorderedLocalParams.index(currentRefVar->GetName());                      // covariances with local scoped
            RooRealVar* currentModVar = dynamic_cast<RooRealVar*>(unorderedLocalParams.at(index)); // fit parameters
            double epsilon = currentRefVar->getError() * 0.01;

            currentModVar->setVal(currentRefVar->getVal() + epsilon);
            double whigh = ReturnCurrentQValue();

            currentModVar->setVal(currentRefVar->getVal() - epsilon);
            double wlow = ReturnCurrentQValue();

            currentModVar->setVal(currentRefVar->getVal());

            derivatives.push_back((whigh - wlow) / (2. * epsilon));
        }

        // Do gaussian error propagation
        const TMatrixDSym cov = rooFitResult->covarianceMatrix();

        double errsq = 0;
        for(int i=0; i<nFreeParams; i++){
            for(int j=0; j<nFreeParams; j++){
                errsq += derivatives.at(i) * cov(i, j) * derivatives.at(j);
            }
        }

        fitResult->weightError = sqrt(errsq);
    }

    // Reset parameters to default values
    for(int i=0; i<nFreeParams; i++){
        RooRealVar* currentRefVar = dynamic_cast<RooRealVar*>(initialParams.at(i));
        int index = unorderedLocalParams.index(currentRefVar->GetName());
        RooRealVar* currentModVar = dynamic_cast<RooRealVar*>(unorderedLocalParams.at(index));
        currentModVar->setVal(currentRefVar->getVal());
    }

    return fitResult;
}



void WibFitFunction::SaveNextFitToFile(std::string fileName){

    saveNextFitToFile = true;
    saveNextFitFileName = fileName;
}



void WibFitFunction::AddData(const PhasespacePoint& phasespacePoint){

    mass->setVal(phasespacePoint.GetMass() - minMass);
    initialWeight->setVal(phasespacePoint.GetInitialWeight());

    if(phasespacePoint.IsMass2Set()){
        mass2->setVal(phasespacePoint.GetMass2() - minMass);
        data->add(RooArgSet(*mass, *mass2, *initialWeight));
        return;
    }

    data->add(RooArgSet(*mass, *initialWeight));
}



double WibFitFunction::GetMinMass() const {

    return minMass;
}



double WibFitFunction::GetMaxMass() const {

    return maxMass;
}
