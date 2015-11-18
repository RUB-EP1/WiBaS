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


#include <iostream>
#include <cmath>
#include <algorithm>

#include "wibas.hh"
#include "fitresult.hh"
#include "phasespace_point.hh"

#include "TTree.h"
#include "TCanvas.h"
#include "TThread.h"
#include "TH2F.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooNumConvPdf.h"
#include "RooProdPdf.h"

const double WiBaS::Pi = 3.1415926;
const bool WiBaS::IS_2PI_CIRCULAR = 1;
const unsigned int WiBaS::FIT_VOIGTIAN = 1;




WiBaS::WiBaS(double pparticleMeanMass, double pparticleWidth, 
	     double pparticleMinMass, double pparticleMaxMass, bool pfit2D) :
   qout(&std::cout),
   fitFunctionType(WiBaS::FIT_VOIGTIAN),
   numNearestNeighbors(200),
   backgroundPolOrder(1),
   particleMeanMass(pparticleMeanMass),
   particleMinMass(pparticleMinMass),
   particleMaxMass(pparticleMaxMass),
   particleWidth(pparticleWidth),
   voigtSigmaMin(particleWidth / 10.),
   voigtSigmaMax(particleWidth * 50.),
   saveNextFitToFile(false),
   fit2D(pfit2D),
   calcErrors(false)
{
   RooMsgService::instance().setSilentMode(true);
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
}



void WiBaS::SetFitFunction(unsigned int pfitFunctionType)
{
   if(fitFunctionType != WiBaS::FIT_VOIGTIAN){
      *qout << "ERROR: fit function currently not supoorted." << std::endl;
      return;
   }

   fitFunctionType = pfitFunctionType;
}



void WiBaS::SetNearestNeighbors(unsigned int pnumNearestNeighbors)
{
   numNearestNeighbors = pnumNearestNeighbors;
}



void WiBaS::SetVoigtianGaussData(double startWidth, double minWidth, double maxWidth){
  voigtSigmaMin = minWidth;
  voigtSigmaMax = maxWidth;
  voigtSigmaStart = startWidth;
}



void WiBaS::SetBackgroundPolOrder(unsigned int order){
  backgroundPolOrder = order;

  if(backgroundPolOrder != 0 && backgroundPolOrder != 1 && backgroundPolOrder != 2){
     *qout << "ERROR: background polinomial order of " << order
	   << " currently not supported." << std::endl;
     exit(0);
  }
}



void WiBaS::RegisterPhasespaceCoord(std::string name, double norm, bool isCircular)
{   
   unsigned short int newID = coordNameMap.size();
   PhasespaceCoord newCoord(newID, norm, isCircular);

   std::pair< std::map< std::string, PhasespaceCoord>::iterator, bool > returnValue;

   returnValue = coordNameMap.insert( std::pair< std::string, PhasespaceCoord >(name, newCoord));
   
   if(returnValue.second == false)
   {
      *qout << "ERROR: element already existing." << std::endl;
   }
}



void WiBaS::AddPhasespacePoint(PhasespacePoint &newPhasespacePoint)
{

   if(!CheckMassInRange(newPhasespacePoint)){
       *qout << "WARNING: Attempt to add a particle outside mass range (m1="
	     << newPhasespacePoint.GetMass() << ", m2="
	     << newPhasespacePoint.GetMass2() << "). " 
 	     << "Rejected." << std::endl;
      return;
   }

   try
   {
      newPhasespacePoint.ArrangeCoordinates(&coordNameMap);
   }
   catch(short int err)
   {
      if(err == PhasespacePoint::ERR_METRIC_MISMATCH)
      {	 
	 *qout << "ERROR: Defined metric has different number of dimensions than the phasespace point."
	       << std::endl;
      }
      else if(err == PhasespacePoint::ERR_UNKNOWN_COORDINATE)
      {
         *qout << "ERROR: Added point contains unknown variables."  << std::endl;
      }
      else *qout << "ERROR: unknown error" << std::endl;
      return;
   }

   PhasespacePoint* copiedPhasespacePoint = new PhasespacePoint(newPhasespacePoint);

   phasespacePointVector.push_back(copiedPhasespacePoint);
}



float WiBaS::CalcPhasespaceDistance(PhasespacePoint* targetPoint, PhasespacePoint* refPoint)
{
   double distance = 0;

   std::map< std::string, PhasespaceCoord>::iterator it;

   for(it=coordNameMap.begin(); it!=coordNameMap.end();++it)
   {
      unsigned short int id = it->second.GetID();
      double norm = it->second.GetNorm();

      if(it->second.GetIsCircular() == true)
      {
	 double distance1 = (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) *
	                    (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) /
	                    (norm * norm);

	 double distance2 = (2*norm - fabs(refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id))) *   
	                    (2*norm - fabs(refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id))) /
	                    (norm * norm);

	 distance += (distance1 < distance2) ? distance1 : distance2; 
      }
      else
      {
	 distance += (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) *
	             (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) /
	             (norm * norm);
      }
   }

   distance = sqrt(distance);
   return distance;
}



void WiBaS::Cleanup()
{
   for(std::vector<PhasespacePoint*>::iterator it = phasespacePointVector.begin();
       it != phasespacePointVector.end(); ++it)
   {
      delete (*it);
   }
   phasespacePointVector.clear();
}



WiBaS::~WiBaS()
{
   Cleanup();
}



bool WiBaS::CalcWeight(PhasespacePoint &refPhasespacePoint)
{
   // check reference point
   try
   {
      refPhasespacePoint.ArrangeCoordinates(&coordNameMap);
   }
   catch(short int err)
   {
      if(err == PhasespacePoint::ERR_METRIC_MISMATCH)
      {	 
	 *qout << "ERROR: Defined metric has different number of dimensions than the phasespace point."
	       << std::endl;
      }
      else if(err == PhasespacePoint::ERR_UNKNOWN_COORDINATE)
      {
         *qout << "ERROR: point of interest contains unknown variables."  << std::endl;
      }
      else *qout << "ERROR: unknown error" << std::endl;
      return false;
   }

   if(!CheckMassInRange(refPhasespacePoint)){
	*qout << "WARNING: Attempt to calculate weight of a particle outside mass range (m1="
	      << refPhasespacePoint.GetMass() << ", m2="
	      << refPhasespacePoint.GetMass2() << "). " << std::endl;
	return false;
   }

   // calculate phasespace distances
   std::vector<FastPointMap> pointMapVector;
   std::vector<PhasespacePoint*>::iterator it;

   for (it=phasespacePointVector.begin(); it!=phasespacePointVector.end(); ++it)
   {   
      FastPointMap newMapEntry((*it), CalcPhasespaceDistance((*it), &refPhasespacePoint));
      pointMapVector.push_back(newMapEntry);
   }

   // sort list
   FastPointMap compHelper(NULL, 0);
   std::sort(pointMapVector.begin(),
   	     pointMapVector.end(), compHelper);

   // Cut vector
   int cutIndex=-1;
   double weightsum=0;
   for(unsigned int i=1; i<pointMapVector.size();i++){
      weightsum += pointMapVector.at(i).phasespacePoint->GetInitialWeight();
      if(weightsum >= numNearestNeighbors){
	 cutIndex = i;
	 break;
      }
   }

   if(cutIndex <= 0){
      *qout << "ERROR: Too few events available for numNearestNeighbors = " << numNearestNeighbors << std::endl;
      return false;
   }
 
   pointMapVector.resize(cutIndex + 1);

   // fill roofit data
   RooRealVar mass("mass","mass", 0, particleMaxMass - particleMinMass); 
   RooRealVar mass2("mass2","mass", 0, particleMaxMass - particleMinMass); 
   RooRealVar initialWeight("initialWeight", "initialWeight", 0);
   RooDataSet* data;

   if(fit2D){
      data = new RooDataSet("data","data", RooArgSet(mass, mass2, initialWeight), RooFit::WeightVar("initialWeight"));
   }
   else{
      data = new RooDataSet("data","data", RooArgSet(mass, initialWeight), RooFit::WeightVar("initialWeight"));
   }

   std::vector<FastPointMap>::iterator it2;
   for(it2=pointMapVector.begin() + 1; it2!=pointMapVector.end(); ++it2) // skip nearest event (=ref event?, TODO: check this!)
   {
      mass.setVal((*it2).phasespacePoint->GetMass() - particleMinMass);
      initialWeight = (*it2).phasespacePoint->GetInitialWeight();

      if(fit2D){
         mass2.setVal((*it2).phasespacePoint->GetMass2() - particleMinMass);
	 data->add(RooArgSet(mass, mass2, initialWeight));

	 mass.setVal((*it2).phasespacePoint->GetMass2() - particleMinMass); // Symmetrization
	 mass2.setVal((*it2).phasespacePoint->GetMass() - particleMinMass);
	 data->add(RooArgSet(mass, mass2, initialWeight));
      }
      else{
	 data->add(RooArgSet(mass, initialWeight));
      }
   } 

   // Do the fit
   FitResult* fitResult;
   if(fitFunctionType == WiBaS::FIT_VOIGTIAN){
      if(fit2D){
	 fitResult = DoVoigtianFit(data, &mass, &mass2,
				   refPhasespacePoint.GetMass() - particleMinMass,
				   refPhasespacePoint.GetMass2() - particleMinMass);
      }
      else{
	 fitResult = DoVoigtianFit(data, &mass, refPhasespacePoint.GetMass() - particleMinMass);
      }
   }
   else{
      *qout << "ERROR: Selected fit function not supported" << std::endl;
      delete data;
      return false;
   }

   // Check fit result
   if((fitResult->rooFitResult == NULL) || (fitResult->rooFitResult->status() != 0))
   {
      *qout << "ERROR: Fit did not converge or returned NULL pointer" << std::endl;
      delete fitResult;
      delete data;
      return false;
   }

   int covQual = fitResult->rooFitResult->covQual();
   if(covQual == 2){
      *qout << "INFO: covariance matrix forced positive-definite" << std::endl;
   }
   else if(covQual == 1){
      *qout << "WARNING: covariance matrix not accurate" << std::endl;
   }
   else if(covQual != 3){
      *qout << "WARNING: covQual = " << covQual << std::endl;
   }
    

   //Get weight at m0
   double Q = fitResult->weight;

   if(Q > 1) 
   {
      *qout << "WARNING: Q > 1. Setting Q = 1." << std::endl;
      Q = 1.0;
   }
   else if(Q < 0)
   {
      *qout << "WARNING: Q < 0. Setting Q = 0." << std::endl;
      Q = 0.0;
   }

   refPhasespacePoint.SetWeight(Q);
   refPhasespacePoint.SetWeightError(fitResult->weightError);

   delete fitResult;
   delete data;

   return true;
}



FitResult* WiBaS::DoVoigtianFit(RooDataSet* data, RooRealVar* mass, double eventMass)
{
   using namespace RooFit;
   double sumOfWeights = data->sumEntries();
   RooRealVar mean("mean","mean / MeV", particleMeanMass - particleMinMass); 
   RooRealVar sigma("sigma","sigma / MeV", voigtSigmaStart, voigtSigmaMin, voigtSigmaMax);
   RooRealVar gamma("gamma","gamma / MeV", particleWidth);
   RooRealVar a1("a1", "a1", 0.1, -100, 100.0);
   RooRealVar a2("a2", "a2", 0.1, -100, 100.0);
   
   RooArgSet bkgArgSet = (backgroundPolOrder == 2) ? RooArgSet(a1,a2) :
                        ((backgroundPolOrder == 1) ? RooArgSet(a1) : RooArgSet());
 
   RooVoigtian voigtFunction("voigt", "signal", *mass, mean, gamma, sigma);
   RooPolynomial polFunction("pol","background", *mass, bkgArgSet);
   RooRealVar sigshare("sigshare","#signal/#total", 0.5, 0, 1);
   RooAddPdf sum("sum","s+b", RooArgList(voigtFunction, polFunction),
                 RooArgList(sigshare));

   RooFitResult* rooFitResult = sum.fitTo(*data, Save(true), 
					  Verbose(false), PrintLevel(-1), PrintEvalErrors(-1));

   FitResult* fitResult = new FitResult;
   fitResult->rooFitResult = rooFitResult;  

   mass->setVal(eventMass);
   RooArgSet invMassArgSet(*mass);
   
   fitResult->weight = QValue(voigtFunction.getVal(&invMassArgSet),
			       polFunction.getVal(&invMassArgSet), sigshare.getVal());


   if(saveNextFitToFile){
      double R = sigshare.getVal();
      saveNextFitToFile = false;
      TCanvas canvas("canvas", "My plots", 0, 0, 550, 500);
      RooPlot *frame = mass->frame();
      data->plotOn(frame, Binning((int)(numNearestNeighbors / 4)), DataError(RooAbsData::SumW2) );
      sum.plotOn(frame, RooFit::Normalization(sumOfWeights, RooAbsReal::NumEvent));
      polFunction.plotOn(frame, RooFit::LineColor(kRed), 
      			 RooFit::Normalization((1-R)*sumOfWeights, RooAbsReal::NumEvent));
      voigtFunction.plotOn(frame, RooFit::LineColor(kGreen), 
      			   RooFit::Normalization(R*sumOfWeights, RooAbsReal::NumEvent));
      frame->Draw();
      canvas.SaveAs(saveFitFileName.c_str());
      delete frame;
   }

   if(!calcErrors){
      return fitResult;
   }

   // Error calculation
   // Get derivatives
   const RooArgList unorderedLocalParams(sigma, a1, a2, sigshare);
   const RooArgList finalParams = rooFitResult->floatParsFinal();
   const int nFreeParams = finalParams.getSize();
   std::vector<double> derivatives;

   for(int i=0; i<nFreeParams; i++){
      RooRealVar* currentRefVar = dynamic_cast<RooRealVar*>(finalParams.at(i));              // Insert better idea to match
      int index = unorderedLocalParams.index(currentRefVar->GetName());                      // covariances with local scoped
      RooRealVar* currentModVar = dynamic_cast<RooRealVar*>(unorderedLocalParams.at(index)); // fit parameters
      double epsilon = currentRefVar->getError() * 0.01;
 
      currentModVar->setVal(currentRefVar->getVal() + epsilon);
      double whigh = QValue(voigtFunction.getVal(&invMassArgSet), polFunction.getVal(&invMassArgSet), sigshare.getVal());

      currentModVar->setVal(currentRefVar->getVal() - epsilon);
      double wlow = QValue(voigtFunction.getVal(&invMassArgSet), polFunction.getVal(&invMassArgSet), sigshare.getVal());

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
   return fitResult;
}  



FitResult* WiBaS::DoVoigtianFit(RooDataSet* data, RooRealVar* mass, RooRealVar* mass2, double eventMass, double eventMass2){
   using namespace RooFit;

   RooRealVar mean("mean","mean / MeV", particleMeanMass - particleMinMass);
   RooRealVar sigma("sigma","sigma / MeV", voigtSigmaStart, voigtSigmaMin, voigtSigmaMax);
   RooRealVar gamma("gamma","gamma / MeV", particleWidth);
   RooRealVar a1("a1", "a1", 0.1, -100, 100.0);
   RooRealVar a2("a2", "a2", 0.1, -100, 100.0);

   RooArgSet bkgArgSet = (backgroundPolOrder == 2) ? RooArgSet(a1,a2) :
                        ((backgroundPolOrder == 1) ? RooArgSet(a1) : RooArgSet());

   RooVoigtian voigtFunction1("signal1", "signal1", *mass, mean, gamma, sigma);
   RooVoigtian voigtFunction2("signal2", "signal2", *mass2, mean, gamma, sigma);

   RooPolynomial polFunction1("background1","background1", *mass, bkgArgSet);
   RooPolynomial polFunction2("background2","background2", *mass2, bkgArgSet);

   RooProdPdf voigtFunction("signal", "signal", voigtFunction1, voigtFunction2);
   RooProdPdf polFunction("background", "background", polFunction1, polFunction2);

   RooRealVar sigshare("sigshare","#signal/#total", 0.5, 0, 1);
   RooAddPdf sum("sum","s+b", RooArgList(voigtFunction, polFunction),
                 RooArgList(sigshare));

   RooFitResult* rooFitResult = sum.fitTo(*data, Save(true),
                                 Verbose(false), PrintLevel(-1), PrintEvalErrors(-1));

   FitResult* fitResult = new FitResult;
   fitResult->rooFitResult = rooFitResult;

   mass->setVal(eventMass);
   mass2->setVal(eventMass2);
   RooArgSet invMassArgSet(*mass, *mass2);

   fitResult->weight = QValue(voigtFunction.getVal(&invMassArgSet),
			       polFunction.getVal(&invMassArgSet), sigshare.getVal());

 
   if(saveNextFitToFile){
      saveNextFitToFile = false;
      TCanvas canvas("canvas", "My plots", 0, 0, 1400, 900);
      canvas.Divide(2,2);

      TH2F* histData = (TH2F*)data->createHistogram("hist", *mass,Binning(9),YVar(*mass2,Binning(9))) ;
      canvas.cd(1);
      histData->Draw("lego");

      TH2F* histVoigt = (TH2F*)voigtFunction.createHistogram("histVoigt", *mass, Binning(100), YVar(*mass2, Binning(100)));
      canvas.cd(2);
      histVoigt->Draw("lego");

      TH2F* histPol = (TH2F*)polFunction.createHistogram("histPol", *mass, Binning(100), YVar(*mass2, Binning(100)));
      canvas.cd(3);
      histPol->Draw("lego");

      TH2F* histFit = (TH2F*)sum.createHistogram("histFit", *mass, Binning(100), YVar(*mass2, Binning(100)));
      canvas.cd(4);
      histFit->Draw("lego");

      canvas.SaveAs(saveFitFileName.c_str());

      delete histData;
      delete histVoigt;
      delete histPol;
      delete histFit;
   }

   if(!calcErrors){
      return fitResult;
   }

   // Error calculation
   // Get derivatives
   const RooArgList unorderedLocalParams(sigma, a1, a2, sigshare);
   const RooArgList finalParams = rooFitResult->floatParsFinal();
   const int nFreeParams = finalParams.getSize();
   std::vector<double> derivatives;

   for(int i=0; i<nFreeParams; i++){
      RooRealVar* currentRefVar = dynamic_cast<RooRealVar*>(finalParams.at(i));
      int index = unorderedLocalParams.index(currentRefVar->GetName());
      RooRealVar* currentModVar = dynamic_cast<RooRealVar*>(unorderedLocalParams.at(index));
      double epsilon = currentRefVar->getError() * 0.01;
 
      currentModVar->setVal(currentRefVar->getVal() + epsilon);
      double whigh = QValue(voigtFunction.getVal(&invMassArgSet), polFunction.getVal(&invMassArgSet), sigshare.getVal());

      currentModVar->setVal(currentRefVar->getVal() - epsilon);
      double wlow = QValue(voigtFunction.getVal(&invMassArgSet), polFunction.getVal(&invMassArgSet), sigshare.getVal());

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
   return fitResult;
}



void WiBaS::SaveNextFitToFile(std::string fileName)
{
   saveNextFitToFile = true;
   saveFitFileName = fileName;
}



bool WiBaS::CheckMassInRange(PhasespacePoint &refPhasespacePoint){

   if(refPhasespacePoint.GetMass() < particleMinMass ||
      refPhasespacePoint.GetMass() > particleMaxMass)
   {
      return false;
   }

   if(fit2D && 
      (refPhasespacePoint.GetMass2() < particleMinMass ||
       refPhasespacePoint.GetMass2() > particleMaxMass))
   {
      return false;
   }
      
   return true;
}



double WiBaS::QValue(double sig, double back, double sigshare){

  double s = sig * sigshare;
  double b = back * (1 - sigshare);

  return s / (s+b);
}



void WiBaS::SetCalcErrors(bool set){
  calcErrors = set;
}


