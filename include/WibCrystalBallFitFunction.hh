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


#ifndef WIBCRYSTALBALLFITFUNCTION_H
#define WIBCRYSTALBALLFITFUNCTION_H

#include "WibFitFunction.hh"

class RooRealVar;
class RooCBShape;
class RooPolynomial;

class WibCrystalBallFitFunction : public WibFitFunction
{
    public:
        WibCrystalBallFitFunction(double particleMeanMass,
                                  double minMass,
                                  double maxMass,
                                  unsigned int backgroundPolOrder,
                                  double sigmaStart,
                                  double sigmaMin,
                                  double sigmaMax,
                                  double alphaStart,
                                  double alphaMin,
                                  double alphaMax,
                                  int pn);

        virtual ~WibCrystalBallFitFunction();
        virtual FitResult* DoFitD(double eventMass, double eventMass2);

    protected:
        virtual double ReturnCurrentQValue();
        virtual void SaveFitToFile(std::string fileName);
        virtual RooArgList GetParamList() const;

    private:
        RooRealVar* mean;
        RooRealVar* sigma;
        RooRealVar* a1;
        RooRealVar* a2;
        RooRealVar* sigshare;
        RooRealVar* alpha;
        RooRealVar* n;

        RooCBShape* cbFunction;
        RooPolynomial* polFunction;
};


#endif
