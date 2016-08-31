{
   gROOT->ProcessLine(".include ../../include");
   gROOT->ProcessLine(".L ../../src/FastPointMap.cc+");
   gROOT->ProcessLine(".L ../../src/PhasespaceCoord.cc+");
   gROOT->ProcessLine(".L ../../src/PhasespacePoint.cc+");
   gROOT->ProcessLine(".L ../../src/PhasespacePointCloud.cc+");
   gROOT->ProcessLine(".L ../../src/WibFitFunction.cc+");
   gROOT->ProcessLine(".L ../../src/WibVoigtFitFunction.cc+");
   gROOT->ProcessLine(".L ../../src/WibasCore.cc+");
   gROOT->ProcessLine(".L testApp.cc+");

   testApp();
}
