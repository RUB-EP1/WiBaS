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

#include "EnergyTest.hh"
#include "PhasespacePoint.hh"



double GaussRand(double width, double mean){

    double u1 = rand() / static_cast<double>(RAND_MAX);
    double u2 = rand() / static_cast<double>(RAND_MAX);
    double x = sqrt(-2. * log(u1)) * cos(2 * 3.14159 * u2) * width + mean;

    return x;
}



int main(){

    // Create two test objects, one to test a bad "fit",
    // a second one to test a good fit
    EnergyTest energyTest1(EnergyTest::DISTANCE_GAUSS);
    EnergyTest energyTest2(EnergyTest::DISTANCE_GAUSS);

    // The phasespace consists of only one variable: x.
    energyTest1.RegisterPhasespaceCoord("x");
    energyTest2.RegisterPhasespaceCoord("x");

    // Create two data sets with size 1000 and two samples from the
    // "Fit result" with size 5000. The events are normally 
    // distributed with width 1 and mean 5, or 5.2 respectively
    // An internal normalization to the distribution variances
    // will be applied automatically
    int nData = 1000;
    int nFit = 5000;
    double mean1 = 5.2;
    double mean2 = 5;
    double width = 1;

    // A random seed is only used for the data generation
    srand((unsigned)time(NULL));

    // Bad fit: the means of the two distributions differ by 0.2
    for(int i=0; i < nData; i++){
        PhasespacePoint newPoint;
        newPoint.SetCoordinate("x", GaussRand(width, mean1));
        energyTest1.AddPhasespacePointData(newPoint);
    }
    for(int i=0; i < nFit; i++){
        PhasespacePoint newPoint;
        newPoint.SetCoordinate("x", GaussRand(width, mean2));
        energyTest1.AddPhasespacePointFit(newPoint);
    }

    // Good fit: both data sets taken from same distribution
    for(int i=0; i < nData; i++){
        PhasespacePoint newPoint;
        newPoint.SetCoordinate("x", GaussRand(width, mean2));
        energyTest2.AddPhasespacePointData(newPoint);
    }
    for(int i=0; i < nFit; i++){
        PhasespacePoint newPoint;
        newPoint.SetCoordinate("x", GaussRand(width, mean2));
        energyTest2.AddPhasespacePointFit(newPoint);
    }

    // Get the reference energies
    double phiFit1 = energyTest1.GetPhi();
    double phiFit2 = energyTest2.GetPhi();
    std::cout << "Phi bad fit = " << phiFit1 << std::endl;
    std::cout << "Phi good fit = " << phiFit2 << std::endl;

    // Get distribution of the test statistics using the
    // resampling technique (multithreaded)   
    // Use only 2*50 points here to limit runtime. The GetResampledPhis
    // method takes a random integer seed > 0 as a optional argument.
    // If no seed is given, one is read vom /dev/urandom
    int resamplingsPerThread = 50;
    int nThreads = 2;
    std::vector<double> resampledPhi1 = energyTest1.GetResampledPhis(resamplingsPerThread, nThreads);
    std::vector<double> resampledPhi2 = energyTest2.GetResampledPhis(resamplingsPerThread, nThreads);

    // Calculate the p-values by comparing the reference energies
    // with the distributions obtained by resampling
    double pValue1 = 0;
    double pValue2 = 0;

    for(auto x : resampledPhi1){
	if(phiFit1 < x)
        pValue1 += 1;
    }
    for(auto x : resampledPhi2){
	if(phiFit2 < x)
	    pValue2 += 1;
    }

    pValue1 /= resampledPhi1.size();
    pValue2 /= resampledPhi2.size();

    // Print the result. Note that occasionally the p-value of the
    // good fit can be poor by chance
    std::cout << "p-value bad fit = " << pValue1 << std::endl;
    std::cout << "p-value good fit = " << pValue2 << std::endl;

    return 0;
}
