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


#ifndef WIBFITFUNCTION_H
#define WIBFITFUNCTION_H

#include <string>
#include "RooArgList.h"
#include "PhasespacePoint.hh"

class RooRealVar;
class RooDataSet;
class RooAddPdf;
class FitResult;

class WibFitFunction
{
public:
  WibFitFunction(double pminMass, double pmaxMass);
  virtual ~WibFitFunction();
  virtual FitResult* DoFitD(double eventMass, double eventMass2)=0;
  void SetCalcErrors(bool set);
  void SaveNextFitToFile(std::string fileName);
  void AddData(const PhasespacePoint& phasespacePoint);
  bool GetCalcError(); 
  double GetMinMass();
  double GetMaxMass();
  FitResult* DoFit(double eventMass, double eventMass2);

protected:
  virtual double ReturnCurrentQValue()=0;
  virtual void SaveFitToFile(std::string fileName)=0;
  virtual RooArgList GetParamList()=0;
  RooRealVar* mass;
  RooRealVar* mass2;
  RooRealVar* initialWeight;
  RooDataSet* data;
  RooAddPdf* totalIntensity;

private:
  bool calcError; 
  bool saveNextFitToFile;
  double minMass;
  double maxMass;
  std::string saveNextFitFileName;
};


#endif
