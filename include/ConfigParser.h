/**
* @file ConfigParser.h
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/
#ifndef _CONFIGPARSER_H_
#define _CONFIGPARSER_H_


#include <TMatrixD.h>
#include <vector>
#include <string>

namespace MDImputation_ns {

class ConfigParser {
  
	public:
  
		/** 
		\brief Class constructor
 		*/
  	ConfigParser();

		/** 
		\brief Class destructor
 		*/
  	virtual ~ConfigParser();
	
		/** 
		\brief Read the config file, parse and set info to be used by other classes
 		*/
		int ReadConfig(std::string filename);


	private:

		/** 
		\brief Print parsed information
 		*/
		void Print();

	public:

		
		//## File info	
		static std::string fConfigFileName;
		static std::string fInputFileName;
		static std::string fOutputFileName;
		//static std::string fConstraintDataFileName;
		//static bool fIsInteractiveRun;
		//static bool fSaveFullInfo;

		//## Algorithm info
		static int fNDim;
		static int fNComponents;
		static int fNIterations;
		//static int fNSimEvents;
	
		
		static bool fUseStoppingCriteria;	
		static double fEpsilon;
		static int fParInitMethod;
		static bool fRandomizeStartPars;
		static bool fRandomizeStartCovariancePars;
		static bool fRandomizeStartMeanPars;
		static bool fFixMeanPars;
		static bool fFixCovariancePars;
		static bool fFixFractionPars;
		static bool fForceDiagonalCovariance;
	
		//## Start user parameters
		static std::vector<TMatrixD> fMeanStartPars;
		static std::vector<TMatrixD> fCovarianceStartPars;
		static std::vector<double> fFractionStartPars;	

		//## Constraint parameters
		static bool fUseConstraints;
		static double fConstraintAlphaScale;
		static double fConstraintAlphaTolerance;
		static bool fUseCovarianceEigenBoundConstraint;
		static std::vector<TMatrixD> fCovarianceEigenMinBound;
		static std::vector<TMatrixD> fCovarianceEigenMaxBound;

		/*
		static bool fUseDeltaLStoppingCriteria;	
		static double fDeltaLEpsilon;
	
		static bool fSetMissingData;
		static double fMissingDataFraction;

		static bool fUseTruncatedEM;
		static bool fUseCensoredEM;
		
	
		static bool fUseKMeansStart;
		static bool fUseRandomStart;
		static bool fUseRandomFromModelStart;	
		static bool fUseRandomMeanStart;
		static bool fUseRandomMeanDiffStart;
		static bool fUseRandomCovarianceStart;
		static bool fUseRandomCovarianceEigenStart;
		static bool fUseRandomDeltaStart;
		static bool fUseRandomNuStart;
		static bool fUseRandomFractionStart;
		static bool fUseMeanImputationStart;

		
		static bool fUseRandomRegenerationAfterStuck;
		
		
		static bool fUseMeanConstraint;
		static bool fUseMeanBoundConstraint;	
		static bool fUseMeanDiffConstraint;
		static bool fUseCovarianceConstraint;
		static bool fUseCovarianceBoundConstraint;
		static bool fUseCovarianceEigenConstraint;
		static bool fUseCovarianceEigenBoundConstraint;
		static bool fUseDeltaBoundConstraint;
		static bool fUseNuBoundConstraint;
		static bool fUseAnnealing;
		static double fAnnealingParStart;	
		static double fAnnealingParStep;
	
		static bool fUseScaleMatrixConstraint;
		static bool fUseScaleMatrixBoundConstraint;
		static bool fUseScaleMatrixEigenConstraint;
		static bool fUseScaleMatrixEigenBoundConstraint;
		static bool fUseLocationConstraint;
		static bool fUseLocationDiffConstraint;
		
		static bool fFixDeltaPar;
		static bool fFixNuPar;
		
		static TMatrixD* fMeanConstraintSign;
		static TMatrixD* fSigmaConstraintSign;
	
		static std::vector<TMatrixD> fMeanStartPar;
		static std::vector<TMatrixD> fCovarianceStartPar;
		static std::vector<TMatrixD> fDeltaStartPar;
		static std::vector<double> fNuStartPar;
		static std::vector<double> fFractionStartPar;

		static std::vector<TMatrixD> fMeanTruePar;
		static std::vector<TMatrixD> fCovarianceTruePar;
		static std::vector<TMatrixD> fDeltaTruePar;
		static std::vector<double> fNuTruePar;
		static std::vector<double> fFractionTruePar;

		static std::vector<TMatrixD> fMeanModelPar;
		static std::vector<TMatrixD> fCovarianceModelPar;
		static std::vector<TMatrixD> fDeltaModelPar;
		static std::vector<double> fNuModelPar;

		
		static bool fUseMeanOffset; 
		static std::vector<TMatrixD> fMeanOffsetPar;

		static bool fUseSigmaOffset; 
		static std::vector<TMatrixD> fSigmaOffsetPar;
	
		static bool fUseDeltaOffset; 
		static std::vector<TMatrixD> fDeltaOffsetPar;


		static std::vector<int> fClusterAGroups;
		static int fNClassificationGroups;
		static std::vector<int> fClusterLnAMin;
		static std::vector<int> fClusterLnAMax;

		static TMatrixD* fMeanParTolerance_min;	
		static TMatrixD* fMeanParTolerance_max;
		static TMatrixD* fMeanDiffParTolerance_min;	
		static TMatrixD* fMeanDiffParTolerance_max;
		static TMatrixD* fSigmaParTolerance_min;
		static TMatrixD* fSigmaParTolerance_max;
		static TMatrixD* fDeltaParTolerance_min;
		static TMatrixD* fDeltaParTolerance_max;

		static TMatrixD* fDataCutMin;	
		static TMatrixD* fDataCutMax;	
	
		static std::vector<double> fDrawMinRange;	
		static std::vector<double> fDrawMaxRange;	
		static std::vector<int> fDrawNBins;	
		*/

	friend class MNMixtureClustering;
		
};//close class

}//close namespace 

#endif
 
