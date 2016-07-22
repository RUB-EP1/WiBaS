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
#include <algorithm>
#include <thread>
#include <sstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

#include "EnergyTest.hh"
#include "PhasespacePoint.hh"




const short EnergyTest::DISTANCE_LOG = 1;
const short EnergyTest::DISTANCE_GAUSS = 2;




EnergyTest::EnergyTest(short distFunc, bool writelog) :
    PhasespacePointCloud(2),
    initialized(false),
    epsilon(1E-6),
    gauss2sigsq(0.04),
    distanceFunc(distFunc)
{
    std::ostringstream filename;

    if(writelog){
        filename << "energygof_class_" << distFunc << ".log";
    }
    else{
        filename << "/dev/null";
    }

    _log.open(filename.str());
    _log << "Energy-Test Goodness of fit evaluator\n\n";

    if(distFunc == DISTANCE_LOG){
        _log << "INFO: using log distance function\n";
    }
    else if(distFunc == DISTANCE_GAUSS){
        _log << "INFO: using gaussian distance function\n";
    }
    else{
        _log << "ERROR: distanceFunction does not exist";
        distanceFunc = DISTANCE_LOG;
    }
}




void EnergyTest::AddPhasespacePointData(PhasespacePoint& newPhasespacePoint){
    PhasespacePointCloud::AddPhasespacePoint(newPhasespacePoint, 1);
}



void EnergyTest::AddPhasespacePointFit(PhasespacePoint& newPhasespacePoint){
    PhasespacePointCloud::AddPhasespacePoint(newPhasespacePoint, 2);
}



double EnergyTest::GetPhi(const std::vector<PhasespacePoint*>& phasespacePointVectorData,
                          const std::vector<PhasespacePoint*>& phasespacePointVectorFit){

    if(!initialized){
        _log << "ERROR: not initialized\n";
        throw;
    }

    double sumOfWeightsData=0;
    double sumOfWeightsFit=0;
    double sumRData=0;
    double sumRDataFit=0;

    int numData = phasespacePointVectorData.size();
    int numFit = phasespacePointVectorFit.size();

    for(int i=0; i<numData;i++){

        sumOfWeightsData += phasespacePointVectorData[i]->GetInitialWeight();

        if(i==numData-1)
            break;

        for(int j=i+1; j<numData;j++){

            if(distanceFunc == EnergyTest::DISTANCE_GAUSS){
                sumRData += phasespacePointVectorData[i]->GetInitialWeight() *
                   phasespacePointVectorData[j]->GetInitialWeight() *
                   RGauss(CalcPhasespaceDistance(phasespacePointVectorData[i],
                                 phasespacePointVectorData[j]));
            }
            else if(distanceFunc == EnergyTest::DISTANCE_LOG){
                sumRData += phasespacePointVectorData[i]->GetInitialWeight() *
                   phasespacePointVectorData[j]->GetInitialWeight() *
                   Rlog(CalcPhasespaceDistance(phasespacePointVectorData[i],
                               phasespacePointVectorData[j]));
            }

        }

    }

    sumRData /= (sumOfWeightsData * sumOfWeightsData);

    for(int i=0; i<numFit;i++){

        sumOfWeightsFit += phasespacePointVectorFit[i]->GetInitialWeight();
        for(int j=0; j<numData;j++){

            if(distanceFunc == EnergyTest::DISTANCE_GAUSS){
                sumRDataFit += phasespacePointVectorFit[i]->GetInitialWeight() *
                   phasespacePointVectorData[j]->GetInitialWeight() *
                   RGauss(CalcPhasespaceDistance(phasespacePointVectorFit[i],
                                 phasespacePointVectorData[j]));
            }
            else if(distanceFunc == EnergyTest::DISTANCE_LOG){
                sumRDataFit += phasespacePointVectorFit[i]->GetInitialWeight() *
                   phasespacePointVectorData[j]->GetInitialWeight() *
                   Rlog(CalcPhasespaceDistance(phasespacePointVectorFit[i],
                               phasespacePointVectorData[j]));
            }
        }
    }

    sumRDataFit /= (sumOfWeightsData * sumOfWeightsFit);

    double phi = sumRData - sumRDataFit;
    _log << "INFO: calculated energy " << phi << "\n";
    _log << "INFO: data weight =  " << sumOfWeightsData << " fit weight = " << sumOfWeightsFit << "\n";

    return phi;
}



double EnergyTest::GetPhi(){
    auto _phasespacePointVectorData = GetPointVector(1);
    auto _phasespacePointVectorFit = GetPointVector(2);
    return GetPhi(_phasespacePointVectorData, _phasespacePointVectorFit);
}



void EnergyTest::CalcNormsFromDistVariances(){

    auto _phasespacePointVectorData = GetPointVector(1);
    auto _phasespacePointVectorFit = GetPointVector(2);
    auto& coordNameMap = GetCoordNameMap();

    double sumOfWeightsData = 0;
    for(auto it = _phasespacePointVectorData.begin(); it != _phasespacePointVectorData.end(); ++it){
        sumOfWeightsData+=(*it)->GetInitialWeight();
    }

    double scaleDataWeights = _phasespacePointVectorData.size() / sumOfWeightsData;
    for(auto it = _phasespacePointVectorData.begin(); it != _phasespacePointVectorData.end(); ++it){
        (*it)->SetInitialWeight((*it)->GetInitialWeight() * scaleDataWeights);
    }


    double sumOfWeightsFit = 0;
    for(auto it = _phasespacePointVectorFit.begin(); it != _phasespacePointVectorFit.end(); ++it){
        sumOfWeightsFit+=(*it)->GetInitialWeight();
    }

    double scaleFitWeights = _phasespacePointVectorFit.size() / sumOfWeightsFit;
    for(auto it = _phasespacePointVectorFit.begin(); it != _phasespacePointVectorFit.end(); ++it){
        (*it)->SetInitialWeight((*it)->GetInitialWeight() * scaleFitWeights);
    }



    for(auto it=coordNameMap.begin(); it!=coordNameMap.end();++it){

        int id = it->second.GetID();
        double meanvalue=0;
        double sumofweights=0;
        for(auto it2 = _phasespacePointVectorFit.begin(); it2 != _phasespacePointVectorFit.end(); ++it2){
            meanvalue += (*it2)->GetCoordValue(id) * (*it2)->GetInitialWeight();
            sumofweights+=(*it2)->GetInitialWeight();
        }
        meanvalue /= sumofweights;

        double sigma=0;
        for(auto it2 = _phasespacePointVectorFit.begin(); it2 != _phasespacePointVectorFit.end(); ++it2){
            sigma += ((*it2)->GetCoordValue(id) - meanvalue) * ((*it2)->GetCoordValue(id) - meanvalue) * (*it2)->GetInitialWeight();
        }
        sigma /= sumofweights;
        sigma = sqrt(sigma);

        if(it->second.GetIsCircular() == true){
            _log << "INFO: Norm of " << (*it).first << " /= 2\n";
            sigma /= 2.;
        }

        it->second.SetNorm(sigma);
        _log << "INFO: Norm of " << (*it).first << " = " << sigma << " Mean = " << meanvalue << "\n";
    }

    double f0max=0;
    for(auto it = _phasespacePointVectorFit.begin(); it != _phasespacePointVectorFit.end(); ++it){
        double weight = (*it)->GetInitialWeight();
        if(weight>f0max)
        f0max = weight;
    }
    epsilon = 1. / (_phasespacePointVectorFit.size() * f0max * 5.);

    _log << "INFO: f0max = " << f0max << " epsilon = " << epsilon << "\n";
    initialized = true;
}



std::vector<double> EnergyTest::GetResampledPhis(long n, short threads){

    if(!initialized){
        _log << "ERROR: not initialized\n";
        throw;
    }

    _log << "INFO: resampling " << n << " times using " << threads << " threads\n";

    ifstream f("/dev/urandom");
    unsigned int seed;
    f.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    srand(seed);
    f.close();

    std::vector<std::thread> theThreads;
    std::vector<std::vector<double> > tPhis;
    tPhis.resize(threads);

    for(int i = 0; i<threads;i++){
        theThreads.push_back(std::thread(&EnergyTest::Threadfunc, this, n, std::ref(tPhis[i])));
    }
    for(auto it = theThreads.begin(); it != theThreads.end(); ++it){
        (*it).join();
    }

    std::vector<double> phis;

    for(auto it = tPhis.begin(); it!=tPhis.end();++it){
        phis.insert(phis.end(), (*it).begin(), (*it).end());
    }

    return phis;
}


void EnergyTest::Threadfunc(long n, std::vector<double>& phis){

    phis.clear();
    std::vector<PhasespacePoint*> phasespacePointVectorData;
    std::vector<PhasespacePoint*> phasespacePointVectorFit;
    std::vector<PhasespacePoint*> phasespacePointVectorTemp;

    auto _phasespacePointVectorData = GetPointVector(1);
    auto _phasespacePointVectorFit = GetPointVector(2);

    double sumofweights = 0;
    for(auto it = _phasespacePointVectorData.begin(); it != _phasespacePointVectorData.end(); ++it){
        sumofweights+=(*it)->GetInitialWeight();
    }

    if(_phasespacePointVectorData.size() + _phasespacePointVectorFit.size() > RAND_MAX){
        std::cerr << "ERROR: RAND_MAX = " << RAND_MAX << std::endl;
    }

    for(int i=0; i<n;i++){
        phasespacePointVectorFit.clear();
        phasespacePointVectorData.clear();
        phasespacePointVectorTemp.clear();

        phasespacePointVectorTemp = _phasespacePointVectorFit;
        phasespacePointVectorTemp.insert(phasespacePointVectorTemp.end(), _phasespacePointVectorData.begin(),
                       _phasespacePointVectorData.end());

        double resamplesumofweights=0;
        while(true){
            long index = rand() % phasespacePointVectorTemp.size();
            resamplesumofweights += phasespacePointVectorTemp.at(index)->GetInitialWeight();

            // use copies for better performance
            PhasespacePoint* newPoint = new PhasespacePoint(*phasespacePointVectorTemp.at(index));
            phasespacePointVectorData.push_back(newPoint);

            phasespacePointVectorTemp.erase(phasespacePointVectorTemp.begin()+index);
            if(resamplesumofweights >= sumofweights){
                break;
            }
        }

        for(unsigned int i=0; i<phasespacePointVectorTemp.size(); i++){
            PhasespacePoint* newPoint = new PhasespacePoint(*phasespacePointVectorTemp.at(i));
            phasespacePointVectorFit.push_back(newPoint);
        }

        phis.push_back(GetPhi(phasespacePointVectorData, phasespacePointVectorFit));

        for(auto it = phasespacePointVectorData.begin(); it != phasespacePointVectorData.end(); ++it){
            delete *it;
        }

        for(auto it = phasespacePointVectorFit.begin(); it != phasespacePointVectorFit.end(); ++it){
            delete *it;
        }
    }
}

