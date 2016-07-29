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

#ifndef PHASESPACEPOINTCLOUD_HH
#define PHASESPACEPOINTCLOUD_HH

#include <string>
#include <map>
#include <vector>

class PhasespacePoint;
class PhasespaceCoord;

class PhasespacePointCloud
{
    public:
        PhasespacePointCloud(int numSubsets=1);
        virtual ~PhasespacePointCloud();
        void RegisterPhasespaceCoord(const std::string& name, double norm=1, bool isCircular=false);

        static const double Pi;
        static const bool IS_2PI_CIRCULAR;

    protected:
        std::ostream* _qout;
        void AddPhasespacePoint(PhasespacePoint& newPhasespacePoint, int subset=1);
        void ArrangePointCoordinates(PhasespacePoint& point);
        float CalcPhasespaceDistance(PhasespacePoint* targetPoint, PhasespacePoint* refPoint);
        std::vector<PhasespacePoint*>& GetPointVector(int subset=1);
        std::map<std::string, PhasespaceCoord>& GetCoordNameMap();

    private:
        std::map< std::string, PhasespaceCoord > _coordNameMap;
        std::vector<std::vector<PhasespacePoint*> > _phasespacePointVectors;

        void Cleanup();
};

#endif // PHASESPACEPOINTCLOUD_HH
