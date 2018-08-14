/**
* @file KMeansClustering.h
* @class KMeansClustering
* @brief KMeansClustering
*
* Perform kmeans clustering
* @author S. Riggi
* @date 30/07/2013
*/


#ifndef _KMEANS_CLUSTERING_h
#define _KMEANS_CLUSTERING_h 1

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TStyle.h>

#include <RInside.h>                    // for the embedded R via RInside


#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>


namespace MDImputation_ns {


class KMeansClustering : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    KMeansClustering();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~KMeansClustering();
	

	public:
		/**
		* \brief Run kmedians clustering
		*/
		int RunKMedians(TMatrixD* dataMatrix,int nComponents,int maxIter=100,int nRandomInit=100);
		/**
		* \brief Run kmedians clustering using an already loaded R table
		*/	
		int RunKMedians(std::string RDataName,int nComponents,int maxIter=100,int nRandomInit=100);


		/**
		* \brief Run clustering
		*/
		int RunKMeans(TMatrixD* dataMatrix,int nComponents,int maxIter=100,int nRandomInit=100);
		/**
		* \brief Run clustering using an already loaded R table
		*/	
		int RunKMeans(std::string RDataName,int nComponents,int maxIter=100,int nRandomInit=100);

	public:
		/**
		* \brief Return cluster sizes 
		*/
		TMatrixD* GetClusterSizes(){return fClusterSizes;}
		/**
		* \brief Return cluster centers 
		*/
		TMatrixD* GetClusterCenters(){return fClusterCenters;}
		/**
		* \brief Return cluster ids 
		*/
		TMatrixD* GetClusterIds(){return fClusterIds;}
		/**
		* \brief Return cluster cov matrix
		*/
		std::vector<TMatrixD*>& GetClusterCovMatrixes(){return fClusterCovMatrixes;}
			
			
	private:
		/**
		* \brief Clear allocated data
		*/
		void ClearData();
		
	private:
		
		TMatrixD* fClusterIds;//list of cluster membership per each data (size=Nx1) 
		TMatrixD* fClusterCenters;//Cluster centers (size=nComponents x Ndim)
		std::vector<TMatrixD*> fClusterCovMatrixes;//Covariance matrix of data per each cluster
		TMatrixD* fClusterSizes;//number of data observations per cluster (nComponents x 1)

		ClassDef(KMeansClustering,1)


};//close class

}//close namespace 

#endif
