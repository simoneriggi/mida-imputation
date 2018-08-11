/**
* @file Util.cc
* @class Util
* @brief Util functions for 
*
* Useful math functions
* @author S. Riggi
* @date 17/09/2012
*/

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


ClassImp(MDImputation_ns::Util)

namespace MDImputation_ns {

RInside Util::fR;

Util::Util(){

}

Util::~Util(){

}

TMatrixD* Util::MakeRandomMissingData(TMatrixD* dataMatrix,double missingDataFraction)
{
	//## Check data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return nullptr;
	}

	//Clone input data
	TMatrixD* dataMatrixWithMiss= (TMatrixD*)dataMatrix->Clone();

	//## Set missing data at random
	long int nDim= dataMatrix->GetNcols();
	long int N= dataMatrix->GetNrows();
	long int NElements= N*nDim;
	double NMissingElements= std::floor(missingDataFraction*NElements);
	
	long int missingCounter= 0;
	std::vector< std::vector<long int> > takenIndex;
	std::vector<long int> alltakenFlag;

	for(long int i=0;i<N;i++){
		takenIndex.push_back( std::vector<long int>() );
		alltakenFlag.push_back(0.);
		for(long int j=0;j<nDim;j++){
			takenIndex[i].push_back(0);
		}//end loop dim
	}//end loop events
	
	while(missingCounter<NMissingElements){

		if(missingCounter%100==0) cout<<"--> "<<missingCounter<<"/"<<NMissingElements<<" events generated..."<<endl;

		long int eventId= static_cast<long int>(gRandom->Uniform(0,N));
		long int variableId= static_cast<long int>(gRandom->Uniform(0,nDim));

		//## Check if this combination has already be taken
		if(takenIndex[eventId][variableId]==0){
			takenIndex[eventId][variableId]= 1;//set as taken

			//## Now check if all patterns have been chosen
			bool isAllTaken= true;
			for(long int j=0;j<nDim;j++){
				if(j!=variableId && takenIndex[eventId][j]== 0) {
					isAllTaken= false;
					break;
				}
			}//end loop dim

			if(!isAllTaken){
				(*dataMatrixWithMiss)(eventId,variableId)= TMath::SignalingNaN();
				missingCounter++;
			}
			else alltakenFlag[eventId]= true; 

		}//close if
	}//end loop while
	
	return dataMatrixWithMiss;

}//close MakeRandomMissingData()


int Util::DumpMatrixToAsciiFile(TMatrixD* dataMatrix,std::string filename)
{
	//Check filename
	if(filename==""){
		cerr<<"WARN: Empty output file specified!"<<endl;
		return -1;
	}

	//Create output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Write to file
	long int N= dataMatrix->GetNrows();
	long int NDim= dataMatrix->GetNcols();
	
	for(long int i=0;i<N;i++){
		for(long int j=0;j<NDim;j++){
			fprintf(fout,"%f  ",(*dataMatrix)(i,j));
		}//end loop dim
		fprintf(fout,"\n");
	}//end loop events

	//Close file
	fclose(fout);

	return 0;

}//close DumpMatrixToAsciiFile()


}//close namespace
