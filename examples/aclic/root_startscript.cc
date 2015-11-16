{
   gROOT->ProcessLine(".L ../../src/phasespace_point.cc+");
   gROOT->ProcessLine(".L ../../src/wibas.cc+");
   gROOT->ProcessLine(".L testApp.cc+");

   testApp();
}
