/**************************************************************
 *                                                            *
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' background suppression                          *
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
 *  This file is part of WiBaS.                               *
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



#ifndef PHASESPACE_POINT_H
#define PHASESPACE_POINT_H

#include <vector>
#include <string>
#include <map>

#include "PhasespaceCoord.hh"




class PhasespacePoint
{
    public:
        PhasespacePoint();
        std::vector<double> coordValueVector;
        std::map< std::string, double > coordValueMap;

        void SetMass(double mass);
        void SetMass2(double mass);
        void SetInitialWeight(double weight);
        void SetCoordinate(std::string name, double value);
        void SetWeight(double weight);
        void SetWeightError(double weightError);
        void ArrangeCoordinates(const std::map< std::string, PhasespaceCoord >& coordNameMap);

        double GetCoordValue(unsigned short int id) const;
        double GetMass() const;
        double GetMass2() const;
        double GetWeight() const;
        double GetWeightError() const;
        double GetInitialWeight() const;

        bool IsMass2Set() const;

        static const short int ERR_METRIC_MISMATCH;
        static const short int ERR_UNKNOWN_COORDINATE;
        static const short int ERR_INDEX_OVERFLOW;

    private:
        double _calculatedEventWeight;
        double _calculatedEventWeightError;
        double _initialWeight;
        double _mass;
        double _mass2;
        double _weight;
        double _weightError;
        bool _mass2Set;
};




#endif
