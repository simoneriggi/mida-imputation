/**
* @file MNMixtureClustering.h
* @class MNMixtureClustering
* @brief MNMixtureClustering
*
* Fit a mixture of multivariate gaussian
* @author S. Riggi
* @date 30/07/2013
*/


#ifndef _MN_MIXTURE_CLUSTERING_h
#define _MN_MIXTURE_CLUSTERING_h 1

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

//Clustering options
struct MNClusteringOptions 
{
	MNClusteringOptions(){
		SetStandardPars();
	};
	~MNClusteringOptions(){};
	
	enum ParInitMethod {
		eKMEANS= 1,
		eUSER= 2,
		eRANDOM= 3
	};
	enum PreImputationMethod {
		eMEAN= 1,
		eLD= 2
	};

	void SetStandardPars()
	{
		//Set standard options
		parInitMethod= eKMEANS;
		randomizeStartPars= false;
		randomizeStartCovariancePars= true;
		randomizeStartMeanPars= true;
		preImputationMethod= eMEAN;
		dataDelimiter= "";
		nComponents= 3;
		nIterations= 100;
		useStoppingCriteria= true;
		epsilon= 5.e-5;
		fixFractionPars= false;
		fixMeanPars= false;
		fixCovariancePars= false;
		forceDiagonalCovariance= false;
		P_start.clear();
		Mu_start.clear();
		Sigma_start.clear();
	}

	//- Options	
	int parInitMethod;
	bool randomizeStartPars;
	bool randomizeStartCovariancePars;
	bool randomizeStartMeanPars;
	int preImputationMethod;
	std::string dataDelimiter;
	int nComponents;
	int nIterations;
	bool useStoppingCriteria;
	double epsilon;
	bool fixFractionPars;
	bool fixMeanPars;
	bool fixCovariancePars;
	bool forceDiagonalCovariance;
	
	//- User mixture pars
	std::vector<double> P_start;
	std::vector<TMatrixD> Mu_start;
	std::vector<TMatrixD> Sigma_start;

};//close MNClusteringOptions struct


class MNMixtureClustering : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MNMixtureClustering();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MNMixtureClustering();
	

	public:
		
		/**
		* \brief Run imputation
		*/
		int RunImputation(std::string filename,MNClusteringOptions options=MNClusteringOptions());

		/**
		* \brief Get imputed data
		*/
		TMatrixD* GetImputedData(){return fDataMatrix_imp;}

	protected:
		/**
		* \brief Run data
		*/
		int ReadData();
		/**
		* \brief Run EM clustering
		*/
		int RunEMClustering();
		/**
		* \brief EM algorithm initialization
		*/
		int RunEM_Init();
		/**
		* \brief Run EM algorithm E-step
		*/
		int RunEM_EStep(double& LL);
		/**
		* \brief Run EM algorithm M-step
		*/
		int RunEM_MStep();
		/**
		* \brief Init mixture parameters to KMeans clustering
		*/
		int InitParsToKMeans();
		/**
		* \brief Init mixture parameters to user provided values
		*/
		int InitParsToUser();
		/**
		* \brief Randomize fit parameters
		*/
		void RandomizePars();
		/**
		* \brief Randomize sigma parameters
		*/
		void RandomizeSigmaPars();
		/**
		* \brief Randomize mean parameters
		*/
		void RandomizeMeanPars();
		/**
		* \brief Print current parameters
		*/
		void PrintPars();
		/**
		* \brief Print current parameters (nicer printout)
		*/
		void PrintPars2();

	private:
		/**
		* \brief Initialize data 
		*/
		int Init();
		/**
		* \brief Delete allocated data
		*/
		void ClearData();
		/**
		* \brief Check covariance matrix for given component
		*/
		int CheckCovariance(int componentId);

	private:

		//- Main variables
		//static int fNComponents;
		static int fNDim;
		static long int fN;

		//- Options	
		static std::string fFileName;
		static MNClusteringOptions fOptions;

		//- Data variables
		std::string fRTableName;
		std::string fRTableName_preImp;
		TMatrixD* fDataMatrix;
		TMatrixD* fDataMatrix_preImp;
		TMatrixD* fDataMatrix_imp;	
		TMatrixD* fMeanData;
		TMatrixD* fCovarianceData;
		TMatrixD* fCovarianceDiagData;

		//static std::vector<TMatrixD> fData;//data (each event of size 1 x fNDim)
		std::vector<TMatrixD> fData;//original data matrix containing missing data (each event of size 1 x fNDim)
		std::vector<TMatrixD> fData_completed;//original data matrix with missing data imputed (each event of size 1 x fNDim)
		std::vector< std::vector<long int> > fObsDataIndexList;
		std::vector< std::vector<long int> > fMissingDataIndexList;
		std::vector<TMatrixD> fData_obs;
		std::vector< std::vector<TMatrixD> > fData_compl;
		std::vector< std::vector<TMatrixD> > fProdSigmaMiss;
		std::vector< std::vector<double> > fTau;

		//- Truncated gaussian correction pars
		std::vector<double> fTruncGaussianNormFactor;
		std::vector<TMatrixD> fMuTruncGaussianCorrection;
		std::vector<TMatrixD> fSigmaTruncGaussianCorrection;

		//- Starting component parameters
		std::vector<double> fP_start;//start mixture weights
		std::vector<TMatrixD> fMu_start;//location for each mixture component of size (1 x fNDim)
		std::vector<TMatrixD> fSigma_start;//covariance for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fSigmaInv_start;
		std::vector<TMatrixD> fSigmaEigen_start;

		//- Boundary pars
		std::vector<TMatrixD> fSigmaEigen_min;
		std::vector<TMatrixD> fSigmaEigen_max;

		//- Fitted component parameters
		static std::vector<double> fP;//mixture weights
		static std::vector<TMatrixD> fMu;//mean for each mixture of size (1 x fNDim)
		static std::vector<TMatrixD> fSigma;//covariance for each mixture of size (fNDim x fNDim)
		static std::vector<TMatrixD> fSigmaInv;
		static std::vector<double> fSigmaDet;
		std::vector<TMatrixD> fSigmaEigen;
		std::vector<TMatrixD> fMuDiff;
		std::vector<TMatrixD> fSigmaDiff;

		//- Safe/fallback component parameters
		std::vector<double> fP_safe;
		std::vector<TMatrixD> fMu_safe;//location for each mixture component of size (1 x fNDim)
		std::vector<TMatrixD> fSigma_safe;//covariance for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fSigmaEigen_safe;//covariance eigenvalues for each mixture of size (fNDim x fNDim)

		//- Fit variables
		double fLogLikelihood;//final LogLikelihood
		std::vector<double> fIterLogLikelihood;
		

	ClassDef(MNMixtureClustering,1)


};//close class

}//close namespace 

#endif
