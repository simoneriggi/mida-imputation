/**
* @file MultipleImputation.cc
* @class MultipleImputation
* @brief MultipleImputation
*
* Perform multiple imputation in data sample
* @author S. Riggi
* @date 18/09/2012
*/

#include <MultipleImputation.h>
#include <Util.h>
#include <DataReader.h>

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
#include <TColor.h>

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

ClassImp(MDImputation_ns::MultipleImputation)

namespace MDImputation_ns {


MultipleImputation::MultipleImputation(){

}

MultipleImputation::~MultipleImputation(){

}


TMatrixD* MultipleImputation::RunImputation(std::string filename,int nRuns,int nRepeatedRuns,std::string delimiter)
{
	//## Check input file
	if(filename==""){
		cerr<<"ERROR: Empty input file specified...exit!"<<endl;
		return nullptr;
	}

	//## Check AMELIA package is present in R and load it
	try{
		Util::fR.parseEvalQ("library(\"Amelia\")");
	}
	catch(...){
		cerr<<"ERROR: Failed to load R package Amelia (hint: check it is installed in your R)!"<<endl;
		return nullptr;
	}

	//## Read data and import to R as a matrix
	cout<<"INFO: Reading data table in R..."<<endl;
	if(DataReader::ReadAsciiInR(filename,delimiter,"dataMI")<0){
		cerr<<"ERROR: Failed to read ascii file and import it in R!"<<endl;
		return nullptr;
	}

	//## Fill missing values with multiple imputation method (MI)
	cout<<"INFO: Fill missing data with MeanImputation (MI) method..."<<endl;
	std::stringstream	ss;
	ss<<"for(i in 1:"<<nRepeatedRuns<<")";
	ss<<"{mi<-amelia(dataMI,m="<<nRuns<<",ts=NULL,cs=NULL);";
	ss<<"if(i==1){dataMatrixMI<-mi$imputations[["<<nRuns<<"]];} else{dataMatrixMI<-dataMatrixMI+mi$imputations[["<<nRuns<<"]];}}";
	std::string RImpCmd= ss.str();
	cout<<"INFO: Running R command: "<<RImpCmd<<endl;
	try{
		//TString RImpCmd= Form("for(i in 1:%d){mi<-amelia(dataMI,m=%d,ts=NULL,cs=NULL);if(i==1){dataMatrixMI<-mi$imputations[[%d]];}else{dataMatrixMI<-dataMatrixMI+mi$imputations[[%d]];}}",fNRepeatedRuns,fNRuns,fNRuns,fNRuns);
		Util::fR.parseEvalQ(std::string(RImpCmd));
		Util::fR.parseEvalQ( std::string(Form("dataMatrixMI <- dataMatrixMI/%d;",nRepeatedRuns)) );
		Util::fR.parseEvalQ("N <- nrow(dataMatrixMI);");
		Util::fR.parseEvalQ("NDim <- ncol(dataMatrixMI);");
	}
	catch(...){
		cerr<<"ERROR: Failures occurred when running multiple imputation inside R!"<<endl;
		return nullptr;
	}

	//## Retrieve results
	TMatrixD* dataMatrixWithImpData= 0;
	try{
		//Rcpp::NumericMatrix dataMatrixMI= Rcpp::as<Rcpp::NumericMatrix>(Util::fR.parseEval("dataMatrixMI"));
		Rcpp::NumericMatrix dataMatrixMI= Util::fR.parseEval("dataMatrixMI");
		long int N= Util::fR.parseEval("N");
		long int NDim= Util::fR.parseEval("NDim");
		dataMatrixWithImpData= new TMatrixD(N,NDim);
		for(long int i=0;i<N;i++){
			for(long int j=0;j<NDim;j++){
				(*dataMatrixWithImpData)(i,j)= dataMatrixMI(i,j);
			}
		}
	}//close try block
	catch(...){
		cerr<<"ERROR: Failed to retrieve data table and relative size with imputed values in R!"<<endl;
		return nullptr;
	}

	return dataMatrixWithImpData;

}//close RunImputation()

TMatrixD* MultipleImputation::RunImputation(TMatrixD* dataMatrix,int nRuns,int nRepeatedRuns)
{
	//## Check data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return nullptr;
	}

	//## WRITE ME!!!
	//...
	//...

	return nullptr;

}//close RunImputation()


}//close namespace 

