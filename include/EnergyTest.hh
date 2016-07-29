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

#ifndef ENERGYTEST_HH
#define ENERGYTEST_HH

#include <fstream>
#include <cmath>

#include "PhasespacePointCloud.hh"
#include "PhasespaceCoord.hh"

class PhasespacePoint;

class EnergyTest : public PhasespacePointCloud
{
    private:
        bool initialized;
        double epsilon;
        double gauss2sigsq;
        short distanceFunc;
        std::ofstream _log;

    public:
        EnergyTest(short distFunc, bool writelog);
        void SetGauss2SigSq(double val){gauss2sigsq = val;}
        double GetPhi();
        double GetPhi(const std::vector<PhasespacePoint*>& phasespacePointVectorData,
                      const std::vector<PhasespacePoint*>& phasespacePointVectorFit);
        std::vector<double> GetResampledPhis(long n, short threads=1);
        void Threadfunc(long n, std::vector<double>& phis);
        void CalcNormsFromDistVariances();
        double Rlog(double distance);
        double RGauss(double distance);
        void AddPhasespacePointData(PhasespacePoint& newPhasespacePoint);
        void AddPhasespacePointFit(PhasespacePoint& newPhasespacePoint);

        static const short DISTANCE_LOG;
        static const short DISTANCE_GAUSS;
};



inline double EnergyTest::Rlog(double distance){

    return -log(distance + epsilon);
}



inline double EnergyTest::RGauss(double distance){

    return exp(-distance*distance / gauss2sigsq);
}


#endif // ENERGYTEST_HH
