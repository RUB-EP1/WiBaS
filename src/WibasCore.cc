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

#include "WibasCore.hh"
#include "FitResult.hh"
#include "PhasespacePoint.hh"
#include "WibFitFunction.hh"

#include "RooMsgService.h"


const double WiBaS::Pi = 3.1415926;
const bool WiBaS::IS_2PI_CIRCULAR = 1;



WiBaS::WiBaS(WibFitFunction& pfitFunction) :
   calcErrors(false),
   numNearestNeighbors(200),
   qout(&std::cout),
   fitFunction(&pfitFunction)
{
   RooMsgService::instance().setSilentMode(true);
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
}



void WiBaS::SetNearestNeighbors(unsigned int pnumNearestNeighbors)
{
   numNearestNeighbors = pnumNearestNeighbors;
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



float WiBaS::CalcPhasespaceDistance(PhasespacePoint* targetPoint, PhasespacePoint* refPoint){
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
   try{
      refPhasespacePoint.ArrangeCoordinates(&coordNameMap);
   }
   catch(short int err){
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

   for (it=phasespacePointVector.begin(); it!=phasespacePointVector.end(); ++it){   
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


   // Fill the fit function with the neighbor data
   std::vector<FastPointMap>::iterator it2;
   for(it2=pointMapVector.begin() + 1; it2!=pointMapVector.end(); ++it2) // skip nearest event (=ref event?, TODO: check this!)
   {
      fitFunction->AddData(*((*it2).phasespacePoint));
   }

			   
   // Do the fit
   FitResult* fitResult = fitFunction->DoFit(refPhasespacePoint.GetMass(),  refPhasespacePoint.GetMass2());


   // Check fit result
   if((fitResult == NULL) ||
      (fitResult->rooFitResult == NULL) || (fitResult->rooFitResult->status() != 0)){
      *qout << "ERROR: Fit did not converge or returned NULL pointer" << std::endl;
      delete fitResult;
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

   if(Q > 1) {
      *qout << "WARNING: Q > 1. Setting Q = 1." << std::endl;
      Q = 1.0;
   }
   else if(Q < 0){
      *qout << "WARNING: Q < 0. Setting Q = 0." << std::endl;
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

  if(refPhasespacePoint.GetMass() < minMass ||
     refPhasespacePoint.GetMass() > maxMass)
  {
     return false;
  }

  if(refPhasespacePoint.IsMass2Set() && 
     (refPhasespacePoint.GetMass2() < minMass ||
      refPhasespacePoint.GetMass2() > maxMass))
  {
     return false;
  }
  
  return true;
}



void WiBaS::SetCalcErrors(bool set){
  calcErrors = set;
  fitFunction->SetCalcErrors(set);
}


