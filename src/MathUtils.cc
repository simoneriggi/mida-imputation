/**
* @file MathUtils.cc
* @class MathUtils
* @brief Mathematical utility functions
*
* Useful math functions
* @author S. Riggi
* @date 17/09/2012
*/


#include <MathUtils.h>
#include <Util.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>
#include <Math/RootFinder.h>

#include <RInside.h>                    // for the embedded R via RInside


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(MDImputation_ns::MathUtils)

namespace MDImputation_ns {


MathUtils::MathUtils(){

}

MathUtils::~MathUtils(){

}


double MathUtils::MGDensityFcn(TMatrixD& y,TMatrixD& mu,TMatrixD& SigmaInv,double SigmaDet)
{
	int dim= SigmaInv.GetNrows();
	//TMatrixD diff(1,dim);
	//diff= TMatrixD(TMatrixD::kTransposed, y-mu);
	TMatrixD diff(dim,1);
	diff= TMatrixD(TMatrixD::kTransposed, y-mu);

	TMatrixD arg(1,1);
	//arg= diff*SigmaInv*(y-mu);
	arg= (y-mu)*SigmaInv*diff;	
	double value= 1./(pow(2*TMath::Pi(),dim/2.)*sqrt(SigmaDet))*exp(-0.5*arg(0,0));

	return value;

}//close MGDensityFcn()



double MathUtils::tauGaus(TMatrixD& y,TMatrixD& mu,TMatrixD& SigmaInv,double SigmaDet,double pi)
{
	double f= MGDensityFcn(y,mu,SigmaInv,SigmaDet);
	double value= pi*f;

	return value;

}//close tauGaus()

}//close namespace 
