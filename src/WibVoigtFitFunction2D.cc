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

#include "WibVoigtFitFunction2D.hh"
#include "FitResult.hh"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooVoigtian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"

#include "TCanvas.h"
#include "TH2F.h"



WibVoigtFitFunction2D::WibVoigtFitFunction2D(double particleMeanMass,
                                             double particleWidth,
                                             double pminMass,
                                             double pmaxMass,
                                             unsigned int backgroundPolOrder,
                                             double voigtSigmaStart,
                                             double voigtSigmaMin,
                                             double voigtSigmaMax) :
    WibFitFunction(pminMass, pmaxMass)
{

    mean = new RooRealVar("mean", "mean", particleMeanMass - pminMass);
    sigma = new RooRealVar("sigma", "sigma", voigtSigmaStart, voigtSigmaMin, voigtSigmaMax);
    gamma = new RooRealVar("gamma", "gamma", particleWidth);
    a1 = new RooRealVar("a1", "a1", 0.1, -100, 100.0);
    a2 = new RooRealVar("a2", "a2", 0.1, -100, 100.);
    sigshare = new RooRealVar("sigshare", "sigshare", 0.5, 0, 1);

    RooArgSet bkgArgSet = (backgroundPolOrder == 2) ?
                           RooArgSet(*a1,*a2) : ((backgroundPolOrder == 1) ? RooArgSet(*a1) : RooArgSet());

    voigtFunction1 = new RooVoigtian("voigt1", "signal1", *mass, *mean, *gamma, *sigma);
    voigtFunction2 = new RooVoigtian("voigt2", "signal2", *mass2, *mean, *gamma, *sigma);
    polFunction1 = new RooPolynomial("background1","background1", *mass, bkgArgSet);
    polFunction2 = new RooPolynomial("background2","background2", *mass2, bkgArgSet);
    voigtFunctionProd = new RooProdPdf("voigt", "signal", *voigtFunction1, *voigtFunction2);
    polFunctionProd = new RooProdPdf("background", "background", *polFunction1, *polFunction2);

    totalIntensity = new RooAddPdf("total","total", RooArgList(*voigtFunctionProd, *polFunctionProd), RooArgList(*sigshare));
    data = new RooDataSet("data","data", RooArgSet(*mass, *mass2, *initialWeight), RooFit::WeightVar("initialWeight"));
}



WibVoigtFitFunction2D::~WibVoigtFitFunction2D(){

    delete mean;
    delete sigma;
    delete gamma;
    delete a1;
    delete a2;
    delete voigtFunctionProd;
    delete polFunctionProd;
    delete voigtFunction1;
    delete voigtFunction2;
    delete polFunction1;
    delete polFunction2;
}



FitResult* WibVoigtFitFunction2D::DoFitD(double eventMass, double eventMass2){
 
    RooFitResult* rooFitResult = totalIntensity->fitTo(*data, RooFit::Save(true),  RooFit::Verbose(false),
                                                       RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1));

    FitResult* fitResult = new FitResult;
    fitResult->rooFitResult = rooFitResult;

    mass->setVal(eventMass);
    fitResult->weight = ReturnCurrentQValue();

    return fitResult;
}



double WibVoigtFitFunction2D::ReturnCurrentQValue(){

    RooArgSet invMassArgSet(*mass);

    double r = sigshare->getVal();
    double s = voigtFunctionProd->getVal(&invMassArgSet) * r;
    double b = polFunctionProd->getVal(&invMassArgSet) * (1 - r);

    return s / (s+b);
}



void WibVoigtFitFunction2D::SaveFitToFile(std::string fileName){

    using namespace RooFit;

    TCanvas canvas("canvas", "My plots", 0, 0, 1400, 900);
    canvas.Divide(2,2);

    TH2F* histData = (TH2F*)data->createHistogram("hist", *mass, Binning(9), YVar(*mass2, Binning(9))) ;
    canvas.cd(1);
    histData->Draw("lego");

    TH2F* histVoigt = (TH2F*)voigtFunctionProd->createHistogram("histVoigt", *mass, Binning(100), YVar(*mass2, Binning(100)));
    canvas.cd(2);
    histVoigt->Draw("lego");

    TH2F* histPol = (TH2F*)polFunctionProd->createHistogram("histPol", *mass, Binning(100), YVar(*mass2, Binning(100)));
    canvas.cd(3);
    histPol->Draw("lego");

    TH2F* histFit = (TH2F*)totalIntensity->createHistogram("histFit", *mass, Binning(100), YVar(*mass2, Binning(100)));
    canvas.cd(4);
    histFit->Draw("lego");

    canvas.SaveAs(fileName.c_str());

    delete histData;
    delete histVoigt;
    delete histPol;
    delete histFit;
}



RooArgList WibVoigtFitFunction2D::GetParamList() const {

    return RooArgList(*sigma, *a1, *a2, *sigshare);
}
