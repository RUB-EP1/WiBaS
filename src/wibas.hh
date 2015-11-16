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

#ifndef WIBAS_H
#define WIBAS_H


#include <string>
#include <map>

#include "phasespace_point.hh"

class FitObj;
class TThread;
class TTree;
class FitObj;
class RooAbsPdf;
class RooRealVar;
class RooFitResult;
class RooDataSet;



class WiBaS
{
 private:
   std::ostream* qout;
   unsigned int fitFunctionType;
   unsigned int numNearestNeighbors;
   unsigned int backgroundPolOrder;
   double particleMeanMass;
   double particleMinMass;
   double particleMaxMass;
   double particleWidth;
   double voigtSigmaStart;
   double voigtSigmaMin;
   double voigtSigmaMax;

   bool saveNextFitToFile;
   bool fit2D;
   std::string saveFitFileName;

   std::map< std::string, PhasespaceCoord > coordNameMap;
   std::vector<PhasespacePoint*> phasespacePointVector;

   float CalcPhasespaceDistance(PhasespacePoint* targetPoint, PhasespacePoint* refPoint);
   FitObj* DoVoigtianFit(RooDataSet* data, RooRealVar* mass, double eventMass);
   FitObj* DoVoigtianFit(RooDataSet* data, RooRealVar* mass, RooRealVar* mass2, double eventMass, double eventMass2);

 public:
   WiBaS(double pparticleMeanMass, double pparticleWidth, 
	 double pparticleMinMass, double pparticleMaxMass, bool pfit2D=false);
   ~WiBaS();
   void SetNearestNeighbors(unsigned int pnumNearestNeighbors);
   void SetFitFunction(unsigned int pfitFunctionType);
   void SetVoigtianGaussData(double startWidth, double minWidth, double maxWidth);
   void SetBackgroundPolOrder(unsigned int order);
   void RegisterPhasespaceCoord(std::string name, double norm, bool isCircular=false);
   void AddPhasespacePoint(PhasespacePoint &newPhasespacePoint);
   void SaveNextFitToFile(std::string fileName);
   void Cleanup();
   double CalcWeight(PhasespacePoint &refPhasespacePoint);
   bool CheckMassInRange(PhasespacePoint &refPhasespacePoint);

   static const double Pi;
   static const bool IS_2PI_CIRCULAR;
   static const unsigned int FIT_VOIGTIAN;
 };



class FitObj
{
 public:
   double weight;
   double weightError;
   RooFitResult* fitResult;

   ~FitObj();   
};




#endif
