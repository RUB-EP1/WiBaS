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


#include <iostream>
#include <stdlib.h>
#include <climits>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "../../src/WibasCore.hh"
#include "../../src/WibVoigtFitFunction.hh"

int testApp(int firstEvent=1, int lastEvent=INT_MAX, bool calcErrors=false)
{
   // Omega mass and width in MeV
   double omegaMass = 782.65;
   double omegaWidth = 8.49;
   double range = 150;


   // Create a voigtian fit function with a polynomial of second order for the background
   WibVoigtFitFunction fitFunction(omegaMass, 
				   omegaWidth, 
				   omegaMass - range, 
				   omegaMass + range,
				   2,    // Order of Background polynomial 
				   10.,  // Voigt's gaussian start width
				   1,    // Voigt's gaussian minimum width
				   30);  // Voigt's gaussian maximum width


   // Initialize the WiBaS object
   WiBaS wibasObj(fitFunction);
   

   // Select the 200 nearest neighbors to determine the background
   wibasObj.SetNearestNeighbors(200);


   // We have three relevant phasespace coordinates: the omega
   // production angle prodTheta and the omega decay angles decTheta and decPhi
   wibasObj.RegisterPhasespaceCoord("prodTheta", 2);  // theta goes vom -1 to 1
   wibasObj.RegisterPhasespaceCoord("decTheta", 2);  // theta goes vom -1 to 1


   // The phi angle goes from -Pi to Pi. However, as this is a circular variable
   // the maximum difference is Pi. This is the way to tell WiBaS:
   wibasObj.RegisterPhasespaceCoord("decPhi", 3.14159, WiBaS::IS_2PI_CIRCULAR);


   // Switch the event error calculation on or off
   wibasObj.SetCalcErrors(calcErrors);


   // Now do some root stuff to load the example data.
   TFile exampleFile("../example.root", "read");
   if(!exampleFile.IsOpen()){
      std::cout << "Data file not found. Execute this program inside the 'standalone' directory" << std::endl;
      exit(0);
   }

   TTree* dataTree = (TTree*)exampleFile.Get("exampletree");
   float prodTheta, decTheta, decPhi, mass;

   dataTree->SetBranchAddress("prodTheta", &prodTheta);
   dataTree->SetBranchAddress("decTheta", &decTheta);
   dataTree->SetBranchAddress("decPhi", &decPhi);
   dataTree->SetBranchAddress("mass", &mass);


   // Get event Range. The optional range is used for
   // parallelization purposes
   firstEvent = std::max(firstEvent, 1);
   lastEvent = std::min(lastEvent, (int)(dataTree->GetEntries()));

   if(firstEvent >= lastEvent){
      std::cout << "ERROR: firstEvent >= lastEvent" << std::endl;
   }
   else{
      std::cout << "Event range: " << firstEvent << " - " << lastEvent << std::endl;
   }


   // Now we need to fill the WiBaS database with ALL available data
   int numTotalEntries = dataTree->GetEntries();
   int numEntriesInRange = (lastEvent - firstEvent + 1);
   for(int i=1; i<=numTotalEntries; i++)
   {
      dataTree->GetEntry(i-1);


      // Create a new point in phasespace and assign coordinates and mass
      PhasespacePoint newPoint;
      newPoint.SetCoordinate("prodTheta", prodTheta);
      newPoint.SetCoordinate("decTheta", decTheta);
      newPoint.SetCoordinate("decPhi", decPhi);
      newPoint.SetMass(mass);
      newPoint.SetMass2(mass);

      // Add it to the WiBaS object
      wibasObj.AddPhasespacePoint(newPoint);
   }


   // In the second run we calculate the weights of all events in the requested
   // range and fill some histograms.
   TH1F* signal = new TH1F("signal", "signal", 100, omegaMass - range, omegaMass + range);
   TH1F* background = new TH1F("background", "background", 100, omegaMass - range, omegaMass + range);
   TH1F* sum = new TH1F("sum", "sum", 100, omegaMass - range, omegaMass + range);
 
   for(int i=firstEvent; i <= lastEvent; i++)
   {
      dataTree->GetEntry(i-1);

      PhasespacePoint newPoint;
      newPoint.SetCoordinate("prodTheta", prodTheta);
      newPoint.SetCoordinate("decTheta", decTheta);
      newPoint.SetCoordinate("decPhi",   decPhi);
      newPoint.SetMass(mass);
      newPoint.SetMass2(mass);

      // Save one example fit
      if(i==10){
	 wibasObj.SaveNextFitToFile("exampleFit.png");
      }


      // Finally: Get the event weight
      wibasObj.CalcWeight(newPoint);
      double Q = newPoint.GetWeight();


      // Status update
      int processedEvents = i-firstEvent+1;
      if(processedEvents % 100 == 0) 
	 std::cout << "Event " << processedEvents << " / " << numEntriesInRange << " (" 
		   << (float)processedEvents/(float)numEntriesInRange * 100.0 << "%)" << std::endl;


      // Fill histograms
      signal->Fill(newPoint.GetMass(), Q);
      background->Fill(newPoint.GetMass(), 1-Q);
      sum->Fill(newPoint.GetMass());
   }


   // Save result histogram
   TCanvas *cResult = new TCanvas("cResult", "cResult", 800, 600);
   sum->SetMinimum(0);
   sum->Draw();
   signal->Draw("same");
   background->Draw("same");
   signal->SetLineColor(kBlue);
   background->SetLineColor(kRed);

   std::ostringstream resultName;
   resultName << "result" << firstEvent << ".png";
   cResult->SaveAs(resultName.str().c_str());

   return 0;
}
