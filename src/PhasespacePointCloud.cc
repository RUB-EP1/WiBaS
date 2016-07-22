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




#include <vector>
#include <iostream>
#include <cmath>

#include "PhasespacePointCloud.hh"
#include "PhasespacePoint.hh"




const double PhasespacePointCloud::Pi = 3.1415926;
const bool PhasespacePointCloud::IS_2PI_CIRCULAR = 1;





PhasespacePointCloud::PhasespacePointCloud(int numSubsets) :
    _qout(&std::cout)
{
    _phasespacePointVectors.resize(numSubsets);
}



PhasespacePointCloud::~PhasespacePointCloud()
{
    Cleanup();
}



void PhasespacePointCloud::Cleanup()
{
    std::vector<std::vector<PhasespacePoint*> >::iterator it1 = _phasespacePointVectors.begin();

    for(it1 = _phasespacePointVectors.begin(); it1 != _phasespacePointVectors.end(); ++it1){

        std::vector<PhasespacePoint*>::iterator it2;

        for(it2 = it1->begin(); it2 != it1->end(); ++it2){
            delete (*it2);
        }
    }

    _phasespacePointVectors.clear();
}



void PhasespacePointCloud::RegisterPhasespaceCoord(const std::string& name, double norm, bool isCircular)
{
    unsigned short int newID = _coordNameMap.size();
    PhasespaceCoord newCoord(newID, norm, isCircular);

    std::pair< std::map< std::string, PhasespaceCoord>::iterator, bool > returnValue;

    returnValue = _coordNameMap.insert( std::pair< std::string, PhasespaceCoord >(name, newCoord));

    if(returnValue.second == false){
        *_qout << "ERROR: element already existing." << std::endl;
    }
}


void PhasespacePointCloud::AddPhasespacePoint(PhasespacePoint &newPhasespacePoint, int subset)
{
    ArrangePointCoordinates(newPhasespacePoint);

    PhasespacePoint* copiedPhasespacePoint = new PhasespacePoint(newPhasespacePoint);

    _phasespacePointVectors.at(subset - 1).push_back(copiedPhasespacePoint);
}



void PhasespacePointCloud::ArrangePointCoordinates(PhasespacePoint& point){

    try{
        point.ArrangeCoordinates(_coordNameMap);
    }

    catch(short int err){
        if(err == PhasespacePoint::ERR_METRIC_MISMATCH){
            *_qout << "ERROR: Defined metric has different number of dimensions than the phasespace point."
                   << std::endl;
        }
        else if(err == PhasespacePoint::ERR_UNKNOWN_COORDINATE)
            *_qout << "ERROR: Added point contains unknown variables."  << std::endl;
        else
            *_qout << "ERROR: unknown error" << std::endl;
        throw err;
    }
}



float PhasespacePointCloud::CalcPhasespaceDistance(PhasespacePoint* targetPoint, PhasespacePoint* refPoint){

    double distance = 0;

    std::map< std::string, PhasespaceCoord>::iterator it;

    for(it=_coordNameMap.begin(); it!=_coordNameMap.end();++it){
        unsigned short int id = it->second.GetID();
        double norm = it->second.GetNorm();

        if(it->second.GetIsCircular() == true){
            double distance1 = (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) *
                               (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) /
                               (norm * norm);

            double distance2 = (2*norm - fabs(refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id))) *
                               (2*norm - fabs(refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id))) /
                               (norm * norm);

            distance += (distance1 < distance2) ? distance1 : distance2;
        }
        else{
            distance += (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) *
                        (refPoint->GetCoordValue(id) - targetPoint->GetCoordValue(id)) /
                        (norm * norm);
        }
    }

    distance = sqrt(distance);
    return distance;
}



std::vector<PhasespacePoint*>& PhasespacePointCloud::GetPointVector(int subset){
    return _phasespacePointVectors.at(subset-1);
}



std::map<std::string, PhasespaceCoord>& PhasespacePointCloud::GetCoordNameMap(){
    return _coordNameMap;
}

