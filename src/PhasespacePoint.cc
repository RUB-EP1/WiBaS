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
#include <math.h>

#include "PhasespacePoint.hh"
#include "WibasCore.hh"


const short int PhasespacePoint::ERR_METRIC_MISMATCH = 1;
const short int PhasespacePoint::ERR_UNKNOWN_COORDINATE = 2;
const short int PhasespacePoint::ERR_INDEX_OVERFLOW = 3;




PhasespacePoint::PhasespacePoint() :
    _calculatedEventWeight(0.),
    _calculatedEventWeightError(0.),
    _initialWeight(1.),
    _mass(0.),
    _mass2(0.),
    _weight(0.),
    _weightError(0.),
    _mass2Set(false)
{
}



void PhasespacePoint::SetCoordinate(std::string name , double value){

    std::pair< std::map<std::string, double>::iterator, bool > returnValue;

    returnValue = coordValueMap.insert( std::pair< std::string, double >(name, value) );

    if(returnValue.second == false){
      // add exception here
    }
}



void PhasespacePoint::SetMass(double mass){

    _mass = mass;
}



void PhasespacePoint::SetMass2(double mass){

    _mass2 = mass;
    _mass2Set = true;
}



bool PhasespacePoint::IsMass2Set() const {

    return _mass2Set;
}



double PhasespacePoint::GetMass() const {

    return _mass;
}



double PhasespacePoint::GetMass2() const {

    return _mass2;
}



void PhasespacePoint::SetWeight(double weight){

    _weight = weight;
}



void PhasespacePoint::SetWeightError(double weightError){

    _weightError = weightError;
}



double PhasespacePoint::GetWeight() const {

    return _weight;
}



double PhasespacePoint::GetWeightError() const{

    return _weightError;
}


double PhasespacePoint::GetInitialWeight() const {

    return _initialWeight;
}



void PhasespacePoint::SetInitialWeight(double weight){

    _initialWeight = weight;
}



void PhasespacePoint::ArrangeCoordinates(const std::map< std::string, PhasespaceCoord >& coordNameMap){

    if(coordNameMap.size() != coordValueMap.size()){
        throw PhasespacePoint::ERR_METRIC_MISMATCH;
        return;
    }

    coordValueVector.clear();
    coordValueVector.resize(coordValueMap.size(), 0);

    for(std::map< std::string, double >::iterator it = coordValueMap.begin(); it != coordValueMap.end(); ++it){

        std::map< std::string, PhasespaceCoord >::const_iterator returnValue = coordNameMap.find(it->first);

        if(returnValue == coordNameMap.end()){
            throw PhasespacePoint::ERR_UNKNOWN_COORDINATE;
            return;
        }

        unsigned short int id = returnValue->second.GetID();
        coordValueVector.at(id) = it->second;
    }
}


double PhasespacePoint::GetCoordValue(unsigned short int id) const {

    if(id >= coordValueVector.size()){
        throw PhasespacePoint::ERR_INDEX_OVERFLOW;
        return 0.0;
    }
    else
        return coordValueVector.at(id);
}



