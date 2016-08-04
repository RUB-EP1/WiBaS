#include <fstream>
#include "Catch-master/single_include/catch.hpp"
#include "WibGaussFitFunction.hh"
#include "RooMsgService.h"
#include "FitResult.hh"



TEST_CASE("WibGaussFitFunction Tests"){

    RooMsgService::instance().setSilentMode(true);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    double mean          = 1000;
    double massmin       = 900;
    double massmax       = 1100;
    double width         = 14;
    double signalShare   = 0.8;
    double polorder      = 1;
    double fitstartwidth = 10;
    double fitminwidth   = 1;
    double fitmaxwidth   = 100;
    int ndata            = 20000;

    WibGaussFitFunction f(mean, massmin, massmax, polorder, 
			  fitstartwidth, fitminwidth, fitmaxwidth);
 
    // Generate Gaussian shape + linear background using Box-Muller transform
    for(int i=0; i<ndata; i++){
        PhasespacePoint newPoint;

        double u1 = rand() / static_cast<double>(RAND_MAX);
        double u2 = rand() / static_cast<double>(RAND_MAX);
        double mass1 = sqrt(-2. * log(u1)) * cos(2 * 3.14159 * u2) * width + mean;
        double mass2 = rand() / static_cast<double>(RAND_MAX) * (massmax - massmin) + massmin;
        double mass = rand() / static_cast<double>(RAND_MAX) < signalShare ? mass1 : mass2;

        newPoint.SetMass(mass);
        f.AddData(newPoint);
    }
    
    std::string pngName = "WibGaussFitFunction.png";
    f.SaveNextFitToFile(pngName);
    FitResult* fitResult = f.DoFit(mean, 0);

    // RooArgList params = f.GetParamList();
    // int index = params.index("mean");
    // REQUIRE((dynamic_cast<RooRealVar*>params.at(i))->getVal() == 14);

    std::ifstream infile(pngName.c_str());
    REQUIRE(infile.good());
    REQUIRE(fitResult->weight == Approx(0.958).epsilon(0.001));

    infile.close();
    delete fitResult;
}
