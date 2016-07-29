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

#include "WibCrystalBallFitFunction.hh"
#include "FitResult.hh"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooAddPdf.h"

#include "TCanvas.h"

WibCrystalBallFitFunction::WibCrystalBallFitFunction(double particleMeanMass,
                                                     double minMass,
                                                     double maxMass,
                                                     unsigned int backgroundPolOrder,
                                                     double sigmaStart,
                                                     double sigmaMin,
                                                     double sigmaMax,
                                                     double alphaStart,
                                                     double alphaMin,
                                                     double alphaMax,
                                                     int pn) :
    WibFitFunction(minMass, maxMass)
{
    mean = new RooRealVar("mean", "mean", particleMeanMass - minMass);
    sigma = new RooRealVar("sigma", "sigma", sigmaStart, sigmaMin, sigmaMax);
    alpha = new RooRealVar("alpha", "alpha", alphaStart, alphaMin, alphaMax);
    n = new RooRealVar("n", "n", pn);
    a1 = new RooRealVar("a1", "a1", 0.1, -100, 100.0);
    a2 = new RooRealVar("a2", "a2", 0.1, -100, 100.);
    sigshare = new RooRealVar("sigshare", "sigshare", 0.5, 0, 1);

    RooArgSet bkgArgSet = (backgroundPolOrder == 2) ? RooArgSet(*a1,*a2) :
    ((backgroundPolOrder == 1) ? RooArgSet(*a1) : RooArgSet());

    cbFunction = new RooCBShape("cbshape", "signal", *mass, *mean, *sigma, *alpha, *n);
    polFunction = new RooPolynomial("background","background", *mass, bkgArgSet);

    totalIntensity = new RooAddPdf("total","total", RooArgList(*cbFunction, *polFunction), RooArgList(*sigshare));
    data = new RooDataSet("data","data", RooArgSet(*mass, *initialWeight), RooFit::WeightVar("initialWeight"));
}



WibCrystalBallFitFunction::~WibCrystalBallFitFunction(){

    delete mean;
    delete sigma;
    delete alpha;
    delete n;
    delete a1;
    delete a2;
    delete cbFunction;
    delete polFunction;
}



FitResult* WibCrystalBallFitFunction::DoFitD(double eventMass, double eventMass2){
 
    RooFitResult* rooFitResult = totalIntensity->fitTo(*data, RooFit::Save(true),  RooFit::Verbose(false),
                                                       RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1));

    FitResult* fitResult = new FitResult;
    fitResult->rooFitResult = rooFitResult;

    mass->setVal(eventMass);
    fitResult->weight = ReturnCurrentQValue();

    return fitResult;
}



double WibCrystalBallFitFunction::ReturnCurrentQValue(){

    RooArgSet invMassArgSet(*mass);

    double r = sigshare->getVal();
    double s = cbFunction->getVal(&invMassArgSet) * r;
    double b = polFunction->getVal(&invMassArgSet) * (1 - r);

    return s / (s+b);
}



void WibCrystalBallFitFunction::SaveFitToFile(std::string fileName){

    double sumOfWeights = data->sumEntries();

    TCanvas canvas("canvas", "My plots", 0, 0, 550, 500);
    RooPlot *frame = mass->frame();

    double r = sigshare->getVal();

    data->plotOn(frame, RooFit::Binning(50), RooFit::DataError(RooAbsData::SumW2) );
    totalIntensity->plotOn(frame, RooFit::Normalization(sumOfWeights, RooAbsReal::NumEvent));
    polFunction->plotOn(frame, RooFit::LineColor(kRed),
                        RooFit::Normalization((1-r)*sumOfWeights, RooAbsReal::NumEvent));
    cbFunction->plotOn(frame, RooFit::LineColor(kGreen),
                       RooFit::Normalization(r*sumOfWeights, RooAbsReal::NumEvent));

    frame->Draw();
    canvas.SaveAs(fileName.c_str());

    delete frame;
}



RooArgList WibCrystalBallFitFunction::GetParamList() const {

    return RooArgList(*sigma, *a1, *a2, *sigshare, *alpha);
}
