#include <MeanImputation.h>
#include <ListwiseDeletion.h>
#include <MultipleImputation.h>
#include <MNMixtureClustering.h>
#include <Util.h>
#include <ConfigParser.h>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TH2.h>
#include <TH3.h>
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
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TPaveText.h>
#include <TVirtualFitter.h>
#include <TObjArray.h>
#include <TMatrixD.h>
#include <TColor.h>
#include <TApplication.h>
#include <TVector3.h>
#include <TView3D.h>
#include <TMarker.h>
#include <TPaletteAxis.h>
#include <TApplication.h>
#include <TGraph2D.h>

#include <RInside.h> // for the embedded R via RInside

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <getopt.h>


using namespace std;
using namespace MDImputation_ns;

enum ImputationMethod {
	eMean= 1,
	eLD= 2,
	eMI= 3,
	eMNClustering= 4,
	eMSNClustering= 5
};


void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input=[INPUT FILENAME] \t Input data file (ascii) with missing data to be imputed "<<endl;
	cout<<"-o, --output=[OUTPUT FILENAME] \t Output data file (ascii) with missing data imputed (default: Output.dat)"<<endl;
	cout<<"-m, --method=[IMPUTATION_METHOD] \t Imputation method (1=Mean,2=Listwise Deletion,3=Multiple imputation,4=Multivariate Normal clustering,5=Multivariate Skew Normal clustering"<<endl;
	cout<<"-r, --nrepeatedruns=[NREPEATED_RUNS] \t Number of repeated runs for MultipleImputation method (default=100)"<<endl;
	cout<<"-R, --nruns=[NRUNS] \t Number of runs for MultipleImputation method (default=5)"<<endl;
	
	cout<<"Clustering Options:"<<endl;
	cout<<"-c, --config=[CONFIG_FILENAME] \t Configuration file with options (NB: override all command line options if provided)"<<endl;
	cout<<"-u, --userpars \t Set starting parameters to user values provided"<<endl;
	cout<<"-d, --ndim=[NDIM] \t Number of variable dimension to be used when setting starting mixture parameters"<<endl;
	cout<<"-k, --ncomponents=[NCOMPONENTS] \t Number of components to be fitted"<<endl;
	cout<<"-M, --means=[CLUSTER_MEANS] \t Starting cluster means "<<endl;
	cout<<"-S, --sigmas=[CLUSTER_COV_MATRIX] \t Starting cluster covariance matrix "<<endl;
	cout<<"-p, --fractions=[CLUSTER_FRACTIONS] \t Starting cluster fractions (must sum unity) "<<endl;
	cout<<"=============================="<<endl;
}

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "output", required_argument, 0, 'o' },
	{ "method", required_argument, 0, 'm' },
	{ "nrepeatedruns", required_argument, 0, 'r' },
	{ "nruns", required_argument, 0, 'R' },
	{ "ncomponents", required_argument, 0, 'k' },
	{ "ndim", required_argument, 0, 'd' },
	{ "means", required_argument, 0, 'M' },
	{ "sigmas", required_argument, 0, 'S' },
	{ "fractions", required_argument, 0, 'p' },
	{ "userpars", no_argument, 0, 'u' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
int imputationMethod= 1;
std::string inputFileName= "";
std::string outputFileName= "Output.dat";
std::string configFileName= "";
long int nRuns= 5;
long int nRepeatedRuns= 100;
int nComponents= 0;
bool setOptionsFromConfigFile= false;
bool useStartingUserPars= false;
int nDim= 0;
std::vector<double> componentWeights_start;
std::vector<double> componentMeans_start;
std::vector<double> componentSigmas_start;

//Functions
int SetMNClusteringOptions(MNClusteringOptions& options);
int SetMNClusteringOptionsFromConfigFile(MNClusteringOptions& options);
TMatrixD* RunMNClustering();

int main(int argc, char **argv)
{
	//====================================================
	//==         PARSE ARGS
	//=====================================================
	//## Check args
	if(argc<2){
		cout<<endl;
		cerr<< "ERROR: Incorrect number of arguments...see program usage!"<<endl;
		Usage(argv[0]);		
		exit(1);
	}

	//## Get command args
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hi:o::m::r::R::k:M::S::p::d::uc::",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
    	case 'i':	
			{
				inputFileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'm':	
			{
				imputationMethod= atoi(optarg);
				break;	
			}
			case 'r':	
			{
				nRepeatedRuns= atol(optarg);
				break;	
			}
			case 'R':	
			{
				nRuns= atol(optarg);
				break;	
			}
			case 'd':	
			{
				nDim= atoi(optarg);
				break;	
			}
			case 'c':	
			{
				configFileName= std::string(optarg);	
				setOptionsFromConfigFile= true;
				break;
			}
			case 'u':	
			{
				useStartingUserPars= true;
				break;	
			}
			case 'k':	
			{
				nComponents= atoi(optarg);
				break;
			}
			case 'p':	
			{
				std::string argStr(optarg);
				Util::ParseStringFields(componentWeights_start,argStr);
				break;
			}
			case 'M':	
			{
				std::string argStr(optarg);
				Util::ParseStringFields(componentMeans_start,argStr);
				break;
			}
			case 'S':	
			{
				std::string argStr(optarg);
				Util::ParseStringFields(componentSigmas_start,argStr);
				break;
			}
    	default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
	
	
	//====================================================
	//==         RUN IMPUTATION
	//=====================================================
	TMatrixD* inputedDataMatrix= 0;
	if(imputationMethod==eMean){
		inputedDataMatrix= MeanImputation::RunImputation(inputFileName,"");
	}
	else if(imputationMethod==eLD){
		inputedDataMatrix= ListwiseDeletion::RunImputation(inputFileName,"");
	}
	else if(imputationMethod==eMI){
		inputedDataMatrix= MultipleImputation::RunImputation(inputFileName,nRuns,nRepeatedRuns);
	}
	else if(imputationMethod==eMNClustering){
		inputedDataMatrix= RunMNClustering();
	}
	else{
		cerr<<"ERROR: Invalid imputation method selected (method="<<imputationMethod<<") ...exit!"<<endl;
		return -1;
	}

	//## Check output
	if(!inputedDataMatrix){
		cerr<<"ERROR: Failed to compute imputed data!"<<endl;
		return -1;
	}

	//## Save data to file
	if(Util::DumpMatrixToAsciiFile(inputedDataMatrix,outputFileName)<0){
		cerr<<"ERROR: Failed to dump imputed data to file!"<<endl;
		return -1;
	}

	return 0;
	
}//close macro


TMatrixD* RunMNClustering()
{
	//## Set options
	MNClusteringOptions options;
	if(setOptionsFromConfigFile){
		if(SetMNClusteringOptionsFromConfigFile(options)<0){
			cerr<<"ERROR: Invalid options passed to MNClustering from config file "<<configFileName<<"!"<<endl;	
			return nullptr;
		}
	}
	else{
		if(SetMNClusteringOptions(options)<0){
			cerr<<"ERROR: Invalid options passed to MNClustering from command line!"<<endl;	
			return nullptr;
		}
	}

	//Run clustering
	MNMixtureClustering* clustering= new MNMixtureClustering;
	if(clustering->RunImputation(inputFileName,options)<0){
		cerr<<"ERROR: MN clustering imputation failed!"<<endl;
		delete clustering;
		clustering= 0;
		return nullptr;
	}
	
	//Retrieve imputed data matrix (copy it because it is destroyed when MNMixtureClustering is destroyed)	
	TMatrixD* inpData= clustering->GetImputedData();
	TMatrixD* inputedDataMatrix= (TMatrixD*)inpData->Clone();

	//Clear data
	if(clustering){
		delete clustering;
		clustering= 0;
	}

	return inputedDataMatrix;

}//close RunMNClustering()

int SetMNClusteringOptionsFromConfigFile(MNClusteringOptions& options)
{
	//Parse config file
	ConfigParser parser;
	if(parser.ReadConfig(configFileName)<0){
		cerr<<"ERROR: Failed to parse config file "<<configFileName<<"!"<<endl;
		return -1;
	}
	
	//Set main options
	inputFileName= ConfigParser::fInputFileName;
	outputFileName= ConfigParser::fOutputFileName;
	options.nComponents= ConfigParser::fNComponents;
	int nDim= ConfigParser::fNDim;
	options.nIterations= ConfigParser::fNIterations;
	options.useStoppingCriteria= ConfigParser::fUseStoppingCriteria;
	options.epsilon= ConfigParser::fEpsilon;

	//Set starting pars options
	options.parInitMethod= ConfigParser::fParInitMethod;
	options.randomizeStartPars= ConfigParser::fRandomizeStartPars;
	options.randomizeStartCovariancePars= ConfigParser::fRandomizeStartCovariancePars;
	options.randomizeStartMeanPars= ConfigParser::fRandomizeStartMeanPars;
	options.P_start= ConfigParser::fFractionStartPars;
	options.Mu_start= ConfigParser::fMeanStartPars;
	options.Sigma_start= ConfigParser::fCovarianceStartPars;

	//Set constraint options
	options.fixFractionPars= ConfigParser::fFixFractionPars;
	options.fixMeanPars= ConfigParser::fFixMeanPars;
	options.fixCovariancePars= ConfigParser::fFixCovariancePars;
	options.forceDiagonalCovariance= ConfigParser::fForceDiagonalCovariance;
	options.useConstraints= ConfigParser::fUseConstraints;
	options.constraintAlphaScale= ConfigParser::fConstraintAlphaScale;
	options.constraintAlphaTolerance= ConfigParser::fConstraintAlphaTolerance;
	options.useCovarianceEigenBoundConstraint= ConfigParser::fUseCovarianceEigenBoundConstraint;
	options.SigmaEigen_min= ConfigParser::fCovarianceEigenMinBound;
	options.SigmaEigen_max= ConfigParser::fCovarianceEigenMaxBound;
	
	return 0;

}//close SetMNClusteringOptionsFromConfigFile()

int SetMNClusteringOptions(MNClusteringOptions& options)
{
	//Check components
	if(nComponents<=0){
		cerr<<"ERROR: Number of components must be >=0!"<<endl;
		return -1;
	}

	//Check if use provided user starting pars
	if(useStartingUserPars){
		options.parInitMethod= MNClusteringOptions::eUSER;

		//Check dim
		if(nDim<=0){
			cerr<<"ERROR: Number of dimensions must be >0!"<<endl;
			return -1;
		}

		//Check fractions weights
		int nFractionPars= static_cast<int>(componentWeights_start.size());
		if(nFractionPars!=nComponents){
			cerr<<"ERROR: Number of components weights given as arg ("<<nFractionPars<<") is different from nComponents ("<<nComponents<<")!"<<endl;
			return -1;
		}
		options.P_start= componentWeights_start;

		//Check component means
		int nMeanPars= static_cast<int>(componentMeans_start.size());
		if(nMeanPars!=nComponents*nDim){
			cerr<<"ERROR: Number of components mean pars given as arg ("<<nMeanPars<<") is different from nComponents x nDim ("<<nComponents*nDim<<")!"<<endl;
			return -1;
		}
		
		//Check sigmas
		int nSigmaPars= static_cast<int>(componentSigmas_start.size());
		if(nSigmaPars != nDim*nDim*nComponents){
			cerr<<"ERROR: Number of components sigma pars given as arg ("<<nSigmaPars<<") is different from nComponents x nDim^2 ("<<nComponents*nDim*nDim<<")!"<<endl;
			return -1;
		}

		//Set mean & sigma pars
		for(int k=0;k<nComponents;k++){
			//TMatrixD meanPars(nDim,1);
			TMatrixD meanPars(1,nDim);
			TMatrixD sigmaPars(nDim,nDim);
			for(int j=0;j<nDim;j++){
				int index= j + k*nDim;	
				//meanPars(j,0)= componentMeans_start[index];		
				meanPars(0,j)= componentMeans_start[index];
				for(int l=0;l<nDim;l++){
					index= k*nDim*nDim + l + j*nDim;
					sigmaPars(j,l)= componentSigmas_start[index];
				}//end loop nDim
			}//end loop nDim
			options.Mu_start.push_back(meanPars);
			options.Sigma_start.push_back(sigmaPars);
			
			cout<<"== COMPONENT "<<k+1<<" =="<<endl;
			cout<<"--> Means"<<endl;
			meanPars.Print();
			cout<<"--> Sigma"<<endl;
			sigmaPars.Print();
			
		}//end loop components

	}//close if

	return 0;

}//close SetMNClusteringOptions()


