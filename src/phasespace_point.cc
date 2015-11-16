/**************************************************************
 *                                                            *            
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' Background Suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2013  Julian Pychy                          *
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

#include "phasespace_point.hh"
#include "wibas.hh"


const short int PhasespacePoint::ERR_METRIC_MISMATCH = 1;
const short int PhasespacePoint::ERR_UNKNOWN_COORDINATE = 2;
const short int PhasespacePoint::ERR_INDEX_OVERFLOW = 3;


PhasespaceCoord::PhasespaceCoord()
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



double PhasespacePoint::GetMass()
{
   return mass;
}



double PhasespacePoint::GetWeight()
{
   return calculatedEventWeight;
}



void PhasespacePoint::SetWeight(double pweight)
{
   calculatedEventWeight = pweight;
}



double PhasespacePoint::GetInitialWeight()
{
   return initialWeight;
}



void PhasespacePoint::SetInitialWeight(double pweight)
{
   initialWeight = pweight;
}



PhasespacePoint::PhasespacePoint()
{
   calculatedEventWeight = 0.0;
   calculatedEventWeightError = 0.0;
   mass = 0.0;
   initialWeight = 1.0;
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



double PhasespacePoint::GetCoordValue(unsigned short int id)
{
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
