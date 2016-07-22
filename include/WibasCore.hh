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

#ifndef WIBASCORE_H
#define WIBASCORE_H


#include <string>
#include <map>
#include <vector>

#include "PhasespacePointCloud.hh"

class FitResult;
class RooAbsPdf;
class RooRealVar;
class RooFitResult;
class RooDataSet;

class WibFitFunction;

class WiBaS : public PhasespacePointCloud
{
    public:
        WiBaS(WibFitFunction& pfitFunction);
        void SetNearestNeighbors(unsigned int pnumNearestNeighbors);
        void AddPhasespacePoint(PhasespacePoint& newPhasespacePoint);
        void SaveNextFitToFile(std::string fileName);
        void SetCalcErrors(bool set=true);
        bool CalcWeight(PhasespacePoint &refPhasespacePoint);

    private:
        bool calcErrors;
        unsigned int numNearestNeighbors;
        WibFitFunction* fitFunction;
        bool CheckMassInRange(PhasespacePoint &refPhasespacePoint) const;
 };






#endif
