/**
* @file DataGenerator.cc
* @class DataGenerator
* @brief Data generator for multivariate normal/skew-normal/skew-t data 
*
* Random generator
* @author S. Riggi
* @date 17/09/2012
*/

#include <DataGenerator.h>
#include <MathUtils.h>
#include <Util.h>
#include <ConfigParser.h>
#include <Logger.h>

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

#include <Math/GSLRndmEngines.h>
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
#include <Math/DistSampler.h>
#include <Math/DistSamplerOptions.h>
#include <Math/Factory.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <climits>

using namespace std;

ClassImp(MDImputation_ns::DataGenerator)

namespace MDImputation_ns {



DataGenerator::DataGenerator(){

	

}//close constructor

DataGenerator::~DataGenerator(){

}

TMatrixD* DataGenerator::GenMNSample (
	int N,
	std::vector<double>& minRange,std::vector<double>& maxRange,
	TMatrixD& mean,TMatrixD& covMatrix,
	bool addMissingData,double missDataFraction
)
{
	//Check mean & covariance matrix
	//...

	//Store list of parameters
	int nDim= covMatrix.GetNrows();
	TMatrixD* dataMatrix= new TMatrixD(N,nDim);

	//Generate data matrix
	int status= GenMNSample(dataMatrix,N,minRange,maxRange,mean,covMatrix,0,addMissingData,missDataFraction);
	if(status<0){
		ERROR_LOG("Failed to generate random sample from multivariate normal!");
		delete dataMatrix;
		dataMatrix= 0;
		return nullptr;
	}

	return dataMatrix;

}//close GenMNSample()


int DataGenerator::GenMNSample (
	TMatrixD* dataMatrix,int N,
	std::vector<double>& minRange,std::vector<double>& maxRange,
	TMatrixD& mean,TMatrixD& covMatrix,
	int dataFillOffset,
	bool addMissingData,double missDataFraction
)
{
	//Check data matrix
	if(!dataMatrix){
		ERROR_LOG("Null ptr to input matrix given!");
		return -1;
	}

	//Check mean & covariance matrix
	if(mean.GetNcols()!=covMatrix.GetNcols()){
		ERROR_LOG("Mean & cov matrix have different sizes!");
		return -1;
	}

	//Check if simmetric matrix
	if(!covMatrix.IsSymmetric() || covMatrix.GetNcols()!=covMatrix.GetNrows()){
		ERROR_LOG("Input covariance matrix is not symmetric and/or square!");
		return -1;
	}

	//Store list of parameters
	int nDim= covMatrix.GetNrows();
	int nPars= 1 + nDim + nDim*(nDim+1)/2;
	
	//Check data matrix size
	if(dataMatrix->GetNrows()<N || dataMatrix->GetNcols()<nDim){
		WARN_LOG("Input data matrix to be fill has size ("<<dataMatrix->GetNrows()<<"x"<<dataMatrix->GetNcols()<<") smaller than data size ("<<N<<"x"<<nDim<<"), resizing it!");
		dataMatrix->ResizeTo(N,nDim);
	}
	
	//Check data range
	std::vector<size_t> dataRangeListSizes {
		minRange.size(), maxRange.size()
	};
	size_t vectSize_min= *std::min_element(dataRangeListSizes.begin(),dataRangeListSizes.end());
	size_t vectSize_max= *std::max_element(dataRangeListSizes.begin(),dataRangeListSizes.end());
	if(vectSize_min==0 || vectSize_min!=vectSize_max || vectSize_min!=nDim){
		ERROR_LOG("Data range vectors have 0 or different sizes or size different from nDim="<<nDim<<"!");
		return -1;
	}

	//Set pars
	double pars[nPars];
	int parCounter= 0;
	pars[parCounter]= 1;
	parCounter++;

	for(int j=0;j<nDim;j++){
		pars[parCounter]= mean(0,j);	
		parCounter++;	
	}//end loop dims

	//Add cov matrix
	for(int j=0;j<nDim;j++){
		for(int l=j;l<nDim;l++){
			pars[parCounter]= covMatrix(j,l);	
			parCounter++;	
		}//end loop dims
	}//end loop dims
		
	cout<<"pars= {";
	for(int i=0;i<nPars;i++) cout<<pars[i]<<",";
	cout<<"}"<<endl;	

	//Create sampling function
	MNPDF mnPDF(nDim);
  TF1* samplingFcn = new TF1("samplingFcn",mnPDF,0,1,nPars);
  samplingFcn->SetParameters(pars);

	//Create and initialize sampler
	ROOT::Math::DistSampler* sampler= ROOT::Math::Factory::CreateDistSampler();
  if(!sampler) {
		WARN_LOG("Default sampler "<<ROOT::Math::DistSamplerOptions::DefaultSampler()<<" is not available trying with Foam..."); 
    ROOT::Math::DistSamplerOptions::SetDefaultSampler("Foam");
  }
  sampler= ROOT::Math::Factory::CreateDistSampler();
  if(!sampler) { 
  	ERROR_LOG("Foam sampler is not available, cannot generate data!");
    return -1;
  }

	sampler->SetFunction(*samplingFcn,nDim);
  sampler->SetRange(minRange.data(),maxRange.data());
  bool ret= sampler->Init();
	if (!ret) { 
  	ERROR_LOG("Error initializing unuran sampler!");
    return -1; 
  }


	//Generate data
	double v[nDim];
	for(int i=0; i<N; ++i) { 
  	sampler->Sample(v);
    for (int j=0; j<nDim;++j){
    	(*dataMatrix)(i+dataFillOffset,j)= v[j];
		}
  }

	//Delete data
	if(samplingFcn){
		delete samplingFcn;
		samplingFcn= 0;
	}
	if(sampler){
		delete sampler;
		sampler= 0;
	}

	//Add missing data?
	if(addMissingData && Util::SetRandomMissingData(dataMatrix,missDataFraction)<0){
		ERROR_LOG("Failed to add missign data to generated data sample!");
		return -1;
	}

	return 0;

}//close GenMNSample()


TMatrixD* DataGenerator::GenMNMixtureSample (
	int N,
	std::vector<double>& minRange,std::vector<double>& maxRange,
	std::vector<double>& weights,std::vector<TMatrixD>& means,std::vector<TMatrixD>& covMatrixes,
	bool addMissingData,double missDataFraction
)
{
	//Check vector sizes
	std::vector<size_t> argListSizes {
		weights.size(), means.size(), covMatrixes.size()
	};
	size_t vectSize_min= *std::min_element(argListSizes.begin(),argListSizes.end());
	size_t vectSize_max= *std::max_element(argListSizes.begin(),argListSizes.end());
	if(vectSize_min==0 || vectSize_min!=vectSize_max){
		ERROR_LOG("Argument vectors (mean/covariance/weights) have 0 or different sizes!");
		return nullptr;
	}
	
	//Store list of parameters
	int nComponents= static_cast<int>(vectSize_min);
	int nDim= covMatrixes[0].GetNrows();
	int nComponentPars= 1 + nDim + nDim*(nDim+1)/2;
	int nPars= nComponentPars*nComponents;
	INFO_LOG("nComponents="<<nComponents<<", nDim="<<nDim<<", nComponentPars="<<nComponentPars<<", nPars="<<nPars);
	
	//Check data range
	std::vector<size_t> dataRangeListSizes {
		minRange.size(), maxRange.size()
	};
	vectSize_min= *std::min_element(dataRangeListSizes.begin(),dataRangeListSizes.end());
	vectSize_max= *std::max_element(dataRangeListSizes.begin(),dataRangeListSizes.end());
	if(vectSize_min==0 || vectSize_min!=vectSize_max || vectSize_min!=nDim){
		ERROR_LOG("Data range vectors have 0 or different sizes or size different from nDim="<<nDim<<"!");
		return nullptr;
	}

	//Set the number of observation to be generated from multinomial 
	ROOT::Math::GSLRngMT* randomEngine = new ROOT::Math::GSLRngMT();
	randomEngine->Initialize();
	std::vector<unsigned int> nEventsPerClass= randomEngine->Multinomial(N,weights);

	//Delete engine
	delete randomEngine;
	randomEngine= 0;

	//Generate data
	int Ntot= 0;
	std::stringstream ss;
	ss<<"nEventsPerClass {";
	for(auto n:nEventsPerClass) {
		Ntot+= n;
		ss<<n<<",";
	}
	ss<<"}";
	INFO_LOG("Ntot="<<Ntot<<", "<<ss.str());
	
	TMatrixD* dataMatrix= new TMatrixD(Ntot,nDim);

	long int offset= 0;
	for(int k=0;k<nComponents;k++){	
		INFO_LOG("Generating sample data for component no. "<<k+1<<"...");
  	if(nEventsPerClass[k]<=0) continue;

		int status= GenMNSample (
			dataMatrix,nEventsPerClass[k],
			minRange,maxRange,
			means[k],covMatrixes[k],
			offset,
			addMissingData,missDataFraction
		);
		offset+= nEventsPerClass[k];

		if(status<0){
			ERROR_LOG("Failed to generate MN sample data for component no. "<<k+1<<"!");
			if(dataMatrix){
				delete dataMatrix;
				dataMatrix= 0;
			}
			return nullptr;
		}//close if
  }//end loop components

	
	return dataMatrix;

}//close GenMNMixtureSample()

/*
TMatrixD* DataGenerator::GenMNMixtureSample (
	int N,
	std::vector<double>& minRange,std::vector<double>& maxRange,
	std::vector<double>& weights,std::vector<TMatrixD>& means,std::vector<TMatrixD>& covMatrixes,
	bool addMissingData,double missDataFraction
)
{
	//Check vector sizes
	std::vector<size_t> argListSizes {
		weights.size(), means.size(), covMatrixes.size()
	};
	size_t vectSize_min= *std::min_element(argListSizes.begin(),argListSizes.end());
	size_t vectSize_max= *std::max_element(argListSizes.begin(),argListSizes.end());
	if(vectSize_min==0 || vectSize_min!=vectSize_max){
		ERROR_LOG("Argument vectors (mean/covariance/weights) have 0 or different sizes!");
		return nullptr;
	}
	
	//Store list of parameters
	int nComponents= static_cast<int>(vectSize_min);
	int nDim= covMatrixes[0].GetNrows();
	int nComponentPars= 1 + nDim + nDim*(nDim+1)/2;
	int nPars= nComponentPars*nComponents;
	INFO_LOG("nComponents="<<nComponents<<", nDim="<<nDim<<", nComponentPars="<<nComponentPars<<", nPars="<<nPars);
	
	//Check data range
	std::vector<size_t> dataRangeListSizes {
		minRange.size(), maxRange.size()
	};
	vectSize_min= *std::min_element(dataRangeListSizes.begin(),dataRangeListSizes.end());
	vectSize_max= *std::max_element(dataRangeListSizes.begin(),dataRangeListSizes.end());
	if(vectSize_min==0 || vectSize_min!=vectSize_max || vectSize_min!=nDim){
		ERROR_LOG("Data range vectors have 0 or different sizes or size different from nDim="<<nDim<<"!");
		return nullptr;
	}

	//Set the number of observation to be generated from multinomial 
	ROOT::Math::GSLRngMT* randomEngine = new GSLRngMT();
	randomEngine->Initialize();
	std::vector<unsigned int> nEventsPerClass= randomEngine->Multinomial(N,weights);
	
	double pars[nPars];
	int parCounter= 0;
	for(int k=0;k<nComponents;k++){
		//Add weights
		pars[parCounter]= weights[k];
		parCounter++;

		//Add means
		for(int j=0;j<nDim;j++){
			pars[parCounter]= means[k](0,j);	
			parCounter++;	
		}//end loop dims

		//Add cov matrix
		for(int j=0;j<nDim;j++){
			for(int l=j;l<nDim;l++){
				pars[parCounter]= covMatrixes[k](j,l);	
				parCounter++;	
			}//end loop dims
		}//end loop dims
		
	}//end loop components

	cout<<"pars= {";
	for(int i=0;i<nPars;i++) cout<<pars[i]<<",";
	cout<<"}"<<endl;
	
	//Create sampling function
	MNMixturePDF mnMixturePDF(nDim,nComponents);
  TF1* samplingFcn = new TF1("samplingFcn",mnMixturePDF,0,1,nPars);
  samplingFcn->SetParameters(pars);

	//Create and initialize sampler
	ROOT::Math::DistSampler* sampler= ROOT::Math::Factory::CreateDistSampler();
  if(!sampler) {
		WARN_LOG("Default sampler "<<ROOT::Math::DistSamplerOptions::DefaultSampler()<<" is not available trying with Foam..."); 
    ROOT::Math::DistSamplerOptions::SetDefaultSampler("Foam");
  }
  sampler= ROOT::Math::Factory::CreateDistSampler();
  if(!sampler) { 
  	ERROR_LOG("Foam sampler is not available, cannot generate data!");
    return nullptr;
  }

	sampler->SetFunction(*samplingFcn,nDim);
  sampler->SetRange(minRange.data(),maxRange.data());
  bool ret= sampler->Init();
	if (!ret) { 
  	ERROR_LOG("Error initializing unuran sampler!");
    return nullptr; 
  }


	//Generate data
	double v[nDim];
	TMatrixD* dataMatrix= new TMatrixD(N,nDim);
	std::vector<double> data(nDim*N);
	for(int i=0; i<N; ++i) { 
  	sampler->Sample(v);
    for (int j=0; j<nDim;++j){
    	(*dataMatrix)(i,j)= v[j];
		}
  }

	//Delete data
	if(samplingFcn){
		delete samplingFcn;
		samplingFcn= 0;
	}
	if(sampler){
		delete sampler;
		sampler= 0;
	}

	//Add missing data?
	if(addMissingData && Util::SetRandomMissingData(dataMatrix,missDataFraction)<0){
		ERROR_LOG("Failed to add missign data to generated data sample!");
		if(dataMatrix){
			delete dataMatrix;
			dataMatrix= 0;
		}	
		return nullptr;
	}

	return dataMatrix;

}//close GenMNMixtureSample()
*/


}//close namespace
