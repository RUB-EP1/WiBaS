/**************************************************************
 *                                                            *            
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' Background Suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2014  Julian Pychy                          *
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
#include <math.h>
#include <algorithm>

#include "wibas.hh"
#include "phasespace_point.hh"

#include "TTree.h"
#include "TCanvas.h"
#include "TThread.h"

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


const double WiBaS::Pi = 3.1415926;
const bool WiBaS::IS_2PI_CIRCULAR = 1;
const unsigned int WiBaS::FIT_VOIGTIAN = 1;
const unsigned int WiBaS::FIT_NOVOSIBIRSK = 2;
const unsigned int WiBaS::FIT_CRYSTALBALL = 3;



WiBaS::WiBaS(double pparticleMeanMass, double pparticleWidth, 
	   double pparticleMinMass, double pparticleMaxMass) :
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
   saveNextFitToFile(false)
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

   if(newPhasespacePoint.GetMass() < particleMinMass ||
      newPhasespacePoint.GetMass() > particleMaxMass)
   {
      *qout << "WARNING: Attempt to add a particle outside mass range (m="
	    << newPhasespacePoint.GetMass() << "). " 
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



void WiBaS::CalcWeight(PhasespacePoint &refPhasespacePoint)
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
      return;
   }

   if(refPhasespacePoint.GetMass() < particleMinMass ||
      refPhasespacePoint.GetMass() > particleMaxMass)
   {
      *qout << "WARNING: calculating weight of particle outside mass range (m="
	    << refPhasespacePoint.GetMass() << "). " << std::endl;
      return;
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
   }
 
   pointMapVector.resize(cutIndex + 1);

   // fill roofit data
   float newMass;
   RooRealVar mass("mass","mass", 0, particleMaxMass - particleMinMass); 
   RooRealVar initialWeight("initialWeight", "initialWeight", 0);
   RooDataSet data("data","data", RooArgSet(mass, initialWeight), RooFit::WeightVar("initialWeight"));
 
   std::vector<FastPointMap>::iterator it2;
   for(it2=pointMapVector.begin() + 1; it2!=pointMapVector.end(); ++it2) // skip nearest event (=ref event?, TODO: check this!)
   {
      newMass = (*it2).phasespacePoint->GetMass() - particleMinMass;
      mass.setVal(newMass);
      data.add(mass, (*it2).phasespacePoint->GetInitialWeight() );
   } 

   // Do the fit
   FitObj* fitResult;
   if(fitFunctionType == WiBaS::FIT_VOIGTIAN){
      fitResult = DoVoigtianFit(&data, &mass, refPhasespacePoint.GetMass() - particleMinMass);
   }
   else{
      *qout << "ERROR: Selected fit function not supported" << std::endl;
      return;
   }

   if((fitResult->fitResult == NULL) || (fitResult->fitResult->status() != 0))
   {
      *qout << "ERROR: Fit did not converge or returned NULL pointer" << std::endl;
      delete fitResult;
      return;
   }

   int covQual = fitResult->fitResult->covQual();
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

   delete fitResult;
}



FitObj* WiBaS::DoVoigtianFit(RooDataSet* data, RooRealVar* mass, double eventMass)
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

   RooFitResult* res = sum.fitTo(*data, Save(true), 
				 Verbose(false), PrintLevel(-1), PrintEvalErrors(-1));

   FitObj* returnFitObj = new FitObj;
   returnFitObj->fitResult = res;  

   mass->setVal(eventMass);
   RooArgSet invMassArgSet(*mass);
   
   double R = sigshare.getVal();
   double Vn = voigtFunction.getVal(&invMassArgSet);
   double s = Vn * R;
   double b = polFunction.getVal(&invMassArgSet) * (1-R);

   returnFitObj->weight = s / (s + b);

   if(saveNextFitToFile)
   {
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

   return returnFitObj;
}  



void WiBaS::SaveNextFitToFile(std::string fileName)
{
   saveNextFitToFile = true;
   saveFitFileName = fileName;
}



FitObj::~FitObj()
{
   delete fitResult;
}
