/**************************************************************
 *                                                            *            
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' Background Suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2014  Julian Pychy                          *
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

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TThread.h"

#include "../../src/wibas.hh"



int testApp()
{   
   #ifdef __CINT__
   std::cout << "This script must be run by ACLiC" << std::endl;
   return 0;
   #endif

   // Omega mass and width in MeV
   double omegaMass = 782.65;
   double omegaWidth = 8.49;


   // We consider a mass window of omega mass +- 150 MeV
   WiBaS wibasObj(omegaMass, omegaWidth, omegaMass - 150, omegaMass + 150);


   // Select the 200 nearest neighbors to determine the background
   wibasObj.SetNearestNeighbors(200);


   // Set fit function to voigtian and set the voigtian's gauss parameters
   // Also set the order of the background polynomial
   wibasObj.SetFitFunction(WiBaS::FIT_VOIGTIAN);
   wibasObj.SetVoigtianGaussData(10., 1., 30); // start width, minimum, maximum
   wibasObj.SetBackgroundPolOrder(1);


   // We have three relevant phasespace coordinates: the omega
   // production angle prodTheta and the omega decay angles decTheta and decPhi
   wibasObj.RegisterPhasespaceCoord("prodTheta", 2);  // theta goes vom -1 to 1
   wibasObj.RegisterPhasespaceCoord("decTheta", 2);  // theta goes vom -1 to 1


   // The phi angle goes from -Pi to Pi. However, as this is a circular variable
   // the maximum difference is Pi. This is the way to tell WiBaS:
   wibasObj.RegisterPhasespaceCoord("decPhi", 3.14159, WiBaS::IS_2PI_CIRCULAR);


   // Now do some root stuff to load the example data.
   TFile exampleFile("../example.root", "read");
   TTree* dataTree = (TTree*)exampleFile.Get("exampletree");
   float prodTheta, decTheta, decPhi, mass;

   dataTree->SetBranchAddress("prodTheta", &prodTheta);
   dataTree->SetBranchAddress("decTheta", &decTheta);
   dataTree->SetBranchAddress("decPhi", &decPhi);
   dataTree->SetBranchAddress("mass", &mass);


   // Now we need to fill the WiBaS database with all available data
   int numEntries = dataTree->GetEntries();

   for(int i=0; i < numEntries; i++)
   {
      dataTree->GetEntry(i);

      // Create a new point in phasespace and assign coordinates and mass
      PhasespacePoint newPoint;
      newPoint.SetCoordinate("prodTheta", prodTheta);
      newPoint.SetCoordinate("decTheta", decTheta);
      newPoint.SetCoordinate("decPhi", decPhi);
      newPoint.SetMass(mass);


      // Add it to the WiBaS object
      wibasObj.AddPhasespacePoint(newPoint);
   }


   // In the second run we calculate the weights of all events
   // and fill some histograms. 
   TH1F* signal = new TH1F("signal", "signal", 100, omegaMass - 150, omegaMass + 150);
   TH1F* background = new TH1F("background", "background", 100, omegaMass - 150, omegaMass + 150);
   TH1F* sum = new TH1F("sum", "sum", 100, omegaMass - 150, omegaMass + 150);

   for(int i=0; i < numEntries; i++)
   {
      dataTree->GetEntry(i);

      PhasespacePoint newPoint;
      newPoint.SetCoordinate("prodTheta", prodTheta);
      newPoint.SetCoordinate("decTheta", decTheta);
      newPoint.SetCoordinate("decPhi", decPhi);
      newPoint.SetMass(mass);

      // Save image for one example fit
      if(i==10){
	 wibasObj.SaveNextFitToFile("exampleFit1.png"); 
      }

      // Finally: Get the event weight
      wibasObj.CalcWeight(newPoint);

      // Status update
      if(i % 100 == 0){
         std::cout << "Event " << i << " / " << numEntries << " (" 
                   << (float)i/(float)numEntries * 100.0 << "%)" << std::endl;
      }

      // Fill histograms
      signal->Fill(newPoint.GetMass(), newPoint.GetWeight());
      background->Fill(newPoint.GetMass(), 1-newPoint.GetWeight());
      sum->Fill(newPoint.GetMass());
   }



   // Save result
   TCanvas *cResult = new TCanvas("cResult", "cResult", 800, 600);
   sum->SetMinimum(0);
   sum->Draw();
   signal->Draw("same");
   background->Draw("same");
   signal->SetLineColor(kBlue);
   background->SetLineColor(kRed);
   cResult->SaveAs("result.png");
   
   return 0;
}


