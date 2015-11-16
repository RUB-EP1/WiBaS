/**************************************************************
 *                                                            *            
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' background suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2012  Julian Pychy                          *
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






class PhasespaceCoord
{
 private:
   unsigned short int id;
   bool isCircular;
   double norm;

 public:
   PhasespaceCoord();
   PhasespaceCoord(unsigned short int pid, double pnorm, bool pisCircular);
   unsigned short int GetID();
   bool GetIsCircular();
   double GetNorm();
};



class PhasespacePoint
{
  private:
   double mass;
   double initialWeight;
   double calculatedEventWeight;
   double calculatedEventWeightError;

  public:
    PhasespacePoint();
    std::vector<double> coordValueVector;
    std::map< std::string, double > coordValueMap;
 
    void SetMass(double pmass);
    void SetWeight(double pweight);
    void SetInitialWeight(double pweight);
    void SetCoordinate(std::string pname, double pvalue);
    void ArrangeCoordinates(std::map< std::string, PhasespaceCoord >* coordNameMap);

    double GetCoordValue(unsigned short int id);
    double GetMass();
    double GetWeight();
    double GetInitialWeight();


    static const short int ERR_METRIC_MISMATCH;
    static const short int ERR_UNKNOWN_COORDINATE;
    static const short int ERR_INDEX_OVERFLOW;
};



class FastPointMap
{
  public:
   FastPointMap();
   FastPointMap(PhasespacePoint* pphasespacePoint, float pdistance);
   PhasespacePoint* phasespacePoint;
   float distance;

   bool operator() (const FastPointMap& i, const FastPointMap& j);
};



#endif