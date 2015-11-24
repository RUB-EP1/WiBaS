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
#include <math.h>

#include "PhasespacePoint.hh"
#include "WibasCore.hh"


const short int PhasespacePoint::ERR_METRIC_MISMATCH = 1;
const short int PhasespacePoint::ERR_UNKNOWN_COORDINATE = 2;
const short int PhasespacePoint::ERR_INDEX_OVERFLOW = 3;



PhasespaceCoord::PhasespaceCoord() :
   id(0), 
   isCircular(false), 
   norm(1)
{

}



PhasespaceCoord::PhasespaceCoord(unsigned short int pid, double pnorm, bool pisCircular) :
   id(pid), 
   isCircular(pisCircular), 
   norm(pnorm)
{

}



unsigned short int PhasespaceCoord::GetID()
{
   return id;
}



bool PhasespaceCoord::GetIsCircular()
{
   return isCircular;
}



double PhasespaceCoord::GetNorm()
{
   return norm;
}



void PhasespacePoint::SetCoordinate(std::string pname , double pvalue)
{

   std::pair< std::map<std::string, double>::iterator, bool > returnValue;

   returnValue = coordValueMap.insert( std::pair< std::string, double >(pname, pvalue) );

   if(returnValue.second == false)
   {
      // add exception here
   }

}



void PhasespacePoint::SetMass(double pmass)
{
   mass = pmass;
}



void PhasespacePoint::SetMass2(double pmass)
{
   mass2 = pmass;
   mass2Set = true;
}



bool PhasespacePoint::IsMass2Set() const {
  return mass2Set;
}



double PhasespacePoint::GetMass() const {
   return mass;
}



double PhasespacePoint::GetMass2() const {
   return mass2;
}



void PhasespacePoint::SetWeight(double pweight){
  weight = pweight;
}



void PhasespacePoint::SetWeightError(double pweightError){
  weightError = pweightError;
}



double PhasespacePoint::GetWeight() const {
  return weight;
}



double PhasespacePoint::GetWeightError() const{
  return weightError;
}


double PhasespacePoint::GetInitialWeight() const {
   return initialWeight;
}



void PhasespacePoint::SetInitialWeight(double pweight)
{
   initialWeight = pweight;
}



PhasespacePoint::PhasespacePoint() :
  calculatedEventWeight(0.),
  calculatedEventWeightError(0.),
  initialWeight(1.),
  mass(0.),
  mass2(0.),
  weight(0.),
  weightError(0.),
  mass2Set(false)
{
}



void PhasespacePoint::ArrangeCoordinates(std::map< std::string, PhasespaceCoord >* coordNameMap)
{
   if(coordNameMap->size() != coordValueMap.size())
   {
      throw PhasespacePoint::ERR_METRIC_MISMATCH;
      return;
   }

   coordValueVector.clear();
   coordValueVector.resize(coordValueMap.size(), 0);
   
   for(std::map< std::string, double >::iterator it = coordValueMap.begin(); 
       it != coordValueMap.end(); ++it)
   {
      std::map< std::string, PhasespaceCoord >::iterator returnValue = coordNameMap->find(it->first);

      if(returnValue == coordNameMap->end())
      {
	 throw PhasespacePoint::ERR_UNKNOWN_COORDINATE;
	 return;
      }

      unsigned short int id = returnValue->second.GetID();
      coordValueVector.at(id) = it->second;
   }

   coordValueMap.clear();
}



double PhasespacePoint::GetCoordValue(unsigned short int id) const {
   if(id >= coordValueVector.size())
   {
      throw PhasespacePoint::ERR_INDEX_OVERFLOW;
      return 0.0;
   }
   else return coordValueVector.at(id);
}



FastPointMap::FastPointMap()
{
   phasespacePoint = NULL;
   distance = 0.0;
}



FastPointMap::FastPointMap(PhasespacePoint* pphasespacePoint, float pdistance) :
   phasespacePoint(pphasespacePoint),
   distance(pdistance)
{
}



bool FastPointMap::operator() (const FastPointMap& i, const FastPointMap& j)
{
   return (i.distance < j.distance);
}
