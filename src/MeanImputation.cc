/**
* @file MeanImputation.cc
* @class MeanImputation
* @brief MeanImputation
*
* Perform mean imputation in data sample
* @author S. Riggi
* @date 18/09/2012
*/

#include <MeanImputation.h>
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

ClassImp(MDImputation_ns::MeanImputation)

namespace MDImputation_ns {

//RInside MeanImputation::fR;

MeanImputation::MeanImputation()
{

}

MeanImputation::~MeanImputation()
{

}
 
int MeanImputation::RunImputationInR(std::string RTableName,std::string RTableName_withoutMiss)
{
	//## Fill missing values with mean imputation method 
	cout<<"INFO: Fill missing data with MeanImputation method..."<<endl;
	
	std::stringstream ss;
	ss<<RTableName_withoutMiss<<" <- apply("<<RTableName<<",2,meanImputation);";
	std::string RCmd= ss.str();
	cout<<"INFO: Running R cmd: "<<RCmd<<endl;

	try{
		Util::fR.parseEvalQ("meanImputation= function(x){x<-as.numeric(as.character(x)); x[is.na(x)]= mean(x, na.rm=TRUE); x;}");
		Util::fR.parseEvalQ(RCmd.c_str());
	}
	catch(...){
		cerr<<"ERROR: Failures occurred when running mean imputation inside R!"<<endl;
		return -1;
	}

	return 0;

}//close RunImputationInR


TMatrixD* MeanImputation::RunImputationFromRTable(std::string RTableName,std::string RTableName_withoutMiss)
{
	//## Fill missing values with mean imputation method 
	cout<<"INFO: Fill missing data with MeanImputation method..."<<endl;
	
	if(RunImputationInR(RTableName,RTableName_withoutMiss)<0){
		cerr<<"ERROR: Failed to impute data of R table "<<RTableName<<" and return imputed R table "<<RTableName_withoutMiss<<"!"<<endl;
		return nullptr;
	}

	//cout<<"DEBUG: Print R table with miss data "<<endl;
	//Util::fR.parseEvalQ(Form("print(%s);",RTableName.c_str()));

	//cout<<"DEBUG: Print R table with miss data imputed"<<endl;
	//Util::fR.parseEvalQ(Form("print(%s);",RTableName_withoutMiss.c_str()));

	//## Retrieve results
	TMatrixD* dataMatrixWithImpData= Util::ConvertRTableToROOTMatrix(RTableName_withoutMiss);
	if(!dataMatrixWithImpData){
		cerr<<"ERROR: Failed to retrieve data table and convert to ROOT matrix!"<<endl;
		return nullptr;
	}

	/*
	std::stringstream ss;
	ss<<"dataMatrix <- apply("<<RTableName<<",2,meanImputation);";
	std::string RCmd= ss.str();
	cout<<"INFO: Running R cmd: "<<RCmd<<endl;

	try{
		Util::fR.parseEvalQ("meanImputation= function(x){x<-as.numeric(as.character(x)); x[is.na(x)]= mean(x, na.rm=TRUE); x;}");
		Util::fR.parseEvalQ(RCmd.c_str());
	}
	catch(...){
		cerr<<"ERROR: Failures occurred when running mean imputation inside R!"<<endl;
		return nullptr;
	}

	//## Retrieve results
	TMatrixD* dataMatrixWithImpData= Util::ConvertRTableToROOTMatrix("dataMatrix");
	if(!dataMatrixWithImpData){
		cerr<<"ERROR: Failed to retrieve data table and convert to ROOT matrix!"<<endl;
		return nullptr;
	}
	*/



	return dataMatrixWithImpData;

}//close RunImputationFromRTable()


TMatrixD* MeanImputation::RunImputation(std::string filename,std::string delimiter)
{
	//## Check input file
	if(filename==""){
		cerr<<"ERROR: Empty input file specified...exit!"<<endl;
		return nullptr;
	}

	//## Read data and import to R as a matrix
	cout<<"INFO: Reading data table in R..."<<endl;
	if(DataReader::ReadAsciiInR(filename,delimiter,"dataMatrix")<0){
		cerr<<"ERROR: Failed to read ascii file and import it in R!"<<endl;
		return nullptr;
	}

	//## Fill missing values with mean imputation method 
	TMatrixD* dataMatrixWithImpData= RunImputationFromRTable("dataMatrix");

	/*
	//## Read data and import to R as a matrix
	cout<<"INFO: Reading data table in R..."<<endl;
	if(DataReader::ReadAsciiInR(filename,delimiter,"dataMI")<0){
		cerr<<"ERROR: Failed to read ascii file and import it in R!"<<endl;
		return nullptr;
	}

	//## Fill missing values with mean imputation method 
	cout<<"INFO: Fill missing data with MeanImputation method..."<<endl;
	try{
		Util::fR.parseEvalQ("meanImputation= function(x){x<-as.numeric(as.character(x)); x[is.na(x)]= mean(x, na.rm=TRUE); x;}");
		Util::fR.parseEvalQ("dataMatrixMI <- apply(dataMI,2,meanImputation);");
		Util::fR.parseEvalQ("N <- nrow(dataMatrixMI);");
		Util::fR.parseEvalQ("NDim <- ncol(dataMatrixMI);");
		//Util::fR.parseEvalQ("print(dataMatrixMI);");
		//Util::fR.parseEvalQ("print(NDim);");
		//Util::fR.parseEvalQ("print(N);");
	}
	catch(...){
		cerr<<"ERROR: Failures occurred when running mean imputation inside R!"<<endl;
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
	*/
	
	return dataMatrixWithImpData;

}//close RunImputation()

TMatrixD* MeanImputation::RunImputation(TMatrixD* dataMatrix)
{
	//## Check data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return nullptr;
	}

	//## Import data matrix in R
	if(Util::ImportMatrixInR(dataMatrix,"dataMatrix")<0){
		cerr<<"ERROR: Failed to import data matrix in R!"<<endl;
		return nullptr;
	}
	
	//## Fill missing values with mean imputation method 
	TMatrixD* dataMatrixWithImpData= RunImputationFromRTable("dataMatrix");


	return dataMatrixWithImpData;

}//close RunImputation()



}//close namespace

