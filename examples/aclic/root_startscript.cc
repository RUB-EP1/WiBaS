{
   gROOT->ProcessLine(".include ../../include");
   gROOT->ProcessLine(".L ../../src/PhasespacePoint.cc+");
   gROOT->ProcessLine(".L ../../src/WibFitFunction.cc+");
   gROOT->ProcessLine(".L ../../src/WibVoigtFitFunction.cc+");
   gROOT->ProcessLine(".L ../../src/WibasCore.cc+");
   gROOT->ProcessLine(".L testApp.cc+");

   testApp();
}
