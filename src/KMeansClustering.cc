/**
* @file KMeansClustering.cc
* @class KMeansClustering
* @brief KMeansClustering
*
* Perform kmeans clustering
* @author S. Riggi
* @date 30/07/2013
*/

#include <KMeansClustering.h>
#include <DataReader.h>
#include <MeanImputation.h>
#include <ListwiseDeletion.h>
#include <MultipleImputation.h>
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

ClassImp(MDImputation_ns::KMeansClustering)

namespace MDImputation_ns {


KMeansClustering::KMeansClustering()
{
	fClusterIds= 0;
	fClusterCenters= 0;
	fClusterSizes= 0;
	fClusterCovMatrixes.clear();
}

KMeansClustering::~KMeansClustering()
{
	//Delete allocated data
	ClearData();

}//close destructor

void KMeansClustering::ClearData()
{
	if(fClusterSizes){
		delete fClusterSizes;
		fClusterSizes= 0;
	}
	if(fClusterIds){
		delete fClusterIds;
		fClusterIds= 0;
	}

	if(fClusterCenters){
		delete fClusterCenters;
		fClusterCenters= 0;
	}

	for(size_t i=0;i<fClusterCovMatrixes.size();i++){
		if(fClusterCovMatrixes[i]){
			delete fClusterCovMatrixes[i];		
			fClusterCovMatrixes[i]= 0;
		}
	}
	fClusterCovMatrixes.clear();

}//close Clear()

int KMeansClustering::RunKMeans(std::string RDataName,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(nComponents<=0){
		cerr<<"ERROR: Invalid numbr of components arg given (must be >0)"<<endl;
		return -1;
	}

	//Clear existing data
	ClearData();
	
	//Define commands to be run in R
	std::stringstream ss;
	ss<<"clusterResults <- kmeans("<<RDataName<<", centers="<<nComponents<<", iter.max="<<maxIter<<", nstart="<<nRandomInit<<");";
	std::string RCmd= ss.str();

	//Run kmeans clustering
	cout<<"INFO: Running kmeans command in R: "<<RCmd<<endl;
	try{
		Util::fR.parseEval(RCmd);
	}
	catch(...){
		cerr<<"ERROR: Kmeans run failed in R!"<<endl;
		return -1;
	}

	//## Retrieve parameters	
	//- Cluster membership: vector of integers (from 1:k) indicating the cluster to which each point is allocated
	cout<<"INFO: Retrieve kmeans clusters ..."<<endl;
	fClusterIds= Util::ConvertRVectToROOTMatrix("clusterResults$cluster");
	if(!fClusterIds){
		cerr<<"ERROR: Failed to retrieve cluster membership vector and convert it to ROOT!"<<endl;
		return -1;
	}

	//- Centers: a matrix of cluster centers
	cout<<"INFO: Retrieve kmeans cluster centers ..."<<endl;
	fClusterCenters= Util::ConvertRTableToROOTMatrix("clusterResults$centers");
	if(!fClusterCenters){
		cerr<<"ERROR: Failed to retrieve cluster center matrix and convert it to ROOT!"<<endl;
		return -1;
	}

	//- Number of data in each cluster
	cout<<"INFO: Retrieve kmeans cluster centers ..."<<endl;
	fClusterSizes= Util::ConvertRVectToROOTMatrix("clusterResults$size");
	if(!fClusterSizes){
		cerr<<"ERROR: Failed to retrieve cluster size vector and convert it to ROOT!"<<endl;
		return -1;
	}

	//## Compute sample covariance matrix for each clusters
	cout<<"INFO: Computing sample covariance of each clusters..."<<endl;
	for(int i=0;i<nComponents;i++){
		//Define cmd
		//cov(data[which(clusterResults$cluster==clusterId),])
		std::string covMatrixName= Form("covMatrix_cluster%d",i+1);
		RCmd= Form("%s <- cov(%s[which(clusterResults$cluster==%d),])",covMatrixName.c_str(),RDataName.c_str(),i+1);
		cout<<"INFO: Running cmd: "<<RCmd<<endl;

		//Exec cmd
		try{
			Util::fR.parseEval(RCmd.c_str());
		}
		catch(...){
			cerr<<"ERROR: Failed to compute covariance matrix for data belonging to cluster no. "<<i+1<<"!"<<endl;
			return -1;				
		}

		//Fill TMatrixD
		TMatrixD* clusterCovMatrix= Util::ConvertRTableToROOTMatrix(covMatrixName);
		if(!clusterCovMatrix){
			cerr<<"ERROR: Failed to retrieve cluster cov matrix and convert it to ROOT!"<<endl;
			return -1;
		}		
		fClusterCovMatrixes.push_back(clusterCovMatrix);
	}//end loop components


	return 0;

}//close RunKMeans()


int KMeansClustering::RunKMeans(TMatrixD* dataMatrix,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return -1;
	}
	if(nComponents<=0){
		cerr<<"ERROR: Invalid numbr of components arg given (must be >0)"<<endl;
		return -1;
	}

	//## Import data matrix in R
	if(Util::ImportMatrixInR(dataMatrix,"dataMatrix")<0){
		cerr<<"ERROR: Failed to import data matrix in R!"<<endl;
		return -1;
	}

	//## Run kmeans on imported data
	return RunKMeans("dataMatrix",nComponents,maxIter,nRandomInit);

}//close RunKMeans()


int KMeansClustering::RunKMedians(std::string RDataName,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(nComponents<=0){
		cerr<<"ERROR: Invalid numbr of components arg given (must be >0)"<<endl;
		return -1;
	}

	//Clear existing data
	ClearData();

	//Load library flexclust
	std::string RLibName= "flexclust";
	if(Util::LoadRLibraries({RLibName})<0){
		cerr<<"ERROR: Failed to load library "<<RLibName<<"!"<<endl;
		return -1;
	}
	
	
	//Define commands to be run in R
	std::stringstream ss;
	ss<<"clusterResults <- kcca("<<RDataName<<",k="<<nComponents<<", family=kccaFamily(\"kmedians\"), save.data=TRUE)";
	std::string RCmd= ss.str();

	//Run kmeans clustering
	cout<<"INFO: Running kmedians command in R: "<<RCmd<<endl;
	try{
		Util::fR.parseEval(RCmd);
	}
	catch(...){
		cerr<<"ERROR: Kmedians run failed in R!"<<endl;
		return -1;
	}

	//## Retrieve parameters	
	//- Cluster membership: vector of integers (from 1:k) indicating the cluster to which each point is allocated
	cout<<"INFO: Retrieve kmedians clusters ..."<<endl;
	fClusterIds= Util::ConvertRVectToROOTMatrix("attributes(clusterResults)$cluster");
	if(!fClusterIds){
		cerr<<"ERROR: Failed to retrieve cluster membership vector and convert it to ROOT!"<<endl;
		return -1;
	}

	//- Centers: a matrix of cluster centers
	cout<<"INFO: Retrieve kmeans cluster centers ..."<<endl;
	fClusterCenters= Util::ConvertRTableToROOTMatrix("attributes(clusterResults)$centers");
	if(!fClusterCenters){
		cerr<<"ERROR: Failed to retrieve cluster center matrix and convert it to ROOT!"<<endl;
		return -1;
	}

	//- Number of data in each cluster
	cout<<"INFO: Retrieve kmeans cluster centers ..."<<endl;
	fClusterSizes= Util::ConvertRVectToROOTMatrix("attributes(clusterResults)$clusinfo$size");
	if(!fClusterSizes){
		cerr<<"ERROR: Failed to retrieve cluster size vector and convert it to ROOT!"<<endl;
		return -1;
	}

	//## Compute sample covariance matrix for each clusters
	cout<<"INFO: Computing sample covariance of each clusters..."<<endl;
	for(int i=0;i<nComponents;i++){
		//Define cmd
		//cov(data[which(clusterResults$cluster==clusterId),])
		std::string covMatrixName= Form("covMatrix_cluster%d",i+1);
		RCmd= Form("%s <- cov(%s[which(attributes(clusterResults)$cluster==%d),])",covMatrixName.c_str(),RDataName.c_str(),i+1);
		cout<<"INFO: Running cmd: "<<RCmd<<endl;

		//Exec cmd
		try{
			Util::fR.parseEval(RCmd.c_str());
		}
		catch(...){
			cerr<<"ERROR: Failed to compute covariance matrix for data belonging to cluster no. "<<i+1<<"!"<<endl;
			return -1;				
		}

		//Fill TMatrixD
		TMatrixD* clusterCovMatrix= Util::ConvertRTableToROOTMatrix(covMatrixName);
		if(!clusterCovMatrix){
			cerr<<"ERROR: Failed to retrieve cluster cov matrix and convert it to ROOT!"<<endl;
			return -1;
		}		
		fClusterCovMatrixes.push_back(clusterCovMatrix);
	}//end loop components


	return 0;

}//close RunKMedians()


int KMeansClustering::RunKMedians(TMatrixD* dataMatrix,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return -1;
	}
	if(nComponents<=0){
		cerr<<"ERROR: Invalid numbr of components arg given (must be >0)"<<endl;
		return -1;
	}

	//## Import data matrix in R
	if(Util::ImportMatrixInR(dataMatrix,"dataMatrix")<0){
		cerr<<"ERROR: Failed to import data matrix in R!"<<endl;
		return -1;
	}

	//## Run kmeans on imported data
	return RunKMedians("dataMatrix",nComponents,maxIter,nRandomInit);

}//close RunKMedians()



}//close namespace

