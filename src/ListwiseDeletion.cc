/**
* @file ListwiseDeletion.cc
* @class ListwiseDeletion
* @brief ListwiseDeletion
*
* Perform list-wise deletion in data sample
* @author S. Riggi
* @date 18/09/2012
*/

#include <ListwiseDeletion.h>
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


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(MDImputation_ns::ListwiseDeletion)

namespace MDImputation_ns {

ListwiseDeletion::ListwiseDeletion(){

}

ListwiseDeletion::~ListwiseDeletion(){

}


TMatrixD* ListwiseDeletion::RunImputation(std::string filename,std::string delimiter)
{
	//## Check input file
	if(filename==""){
		cerr<<"ERROR: Empty input file specified...exit!"<<endl;
		return nullptr;
	}

	//## Read data matrix from file
	TMatrixD* dataMatrix= DataReader::ReadAscii(filename,delimiter);
	if(!dataMatrix){
		cerr<<"ERROR: Failed to read data from file "<<filename<<"!"<<endl;
		return nullptr;
	}

	//## Delete NAN lines from matrix
	TMatrixD* dataMatrixWithoutMissData= RunImputation(dataMatrix);

	if(dataMatrix){
		delete dataMatrix;
		dataMatrix= 0;
	}

	if(!dataMatrixWithoutMissData){
		cerr<<"ERROR: Failed to compute data without missing values!"<<endl;
		return nullptr;
	}

	return dataMatrixWithoutMissData;

}//close RunImputation()

TMatrixD* ListwiseDeletion::RunImputation(TMatrixD* dataMatrix)
{
	//## Check data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return nullptr;
	}

	//## Loop over matrix and select lines without missing data
	long int N= dataMatrix->GetNrows();
	long int NDim= dataMatrix->GetNcols();
	std::vector<long int> selDataIndexes;	

	for(long int i=0;i<N;i++){
		bool hasMissingData= false;
		for(long int j=0;j<NDim;j++){
			double w= (*dataMatrix)(i,j);
			if(TMath::IsNaN(w) || fabs(w)==TMath::Infinity()){
				hasMissingData= true;
				break;
			}
		}//end loop dims
		
		if(!hasMissingData) selDataIndexes.push_back(i);

	}//end loop entries
	
	cout<<"INFO: Selected "<<selDataIndexes.size()<<"/"<<N<<" data entries..."<<endl;

	TMatrixD* dataMatrixWithoutMissData= new TMatrixD(selDataIndexes.size(),NDim);
	
	for(size_t i=0;i<selDataIndexes.size();i++){
		long int index= selDataIndexes[i];
		for(long int j=0;j<NDim;j++){
			double w= (*dataMatrix)(index,j);
			(*dataMatrixWithoutMissData)(i,j)= w;
		}//end loop dims
	}//end loop 
	
	return dataMatrixWithoutMissData;

}//close RunImputation()


}//close namespace

