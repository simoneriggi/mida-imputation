/**
* @file MNMixtureClustering.cc
* @class MNMixtureClustering
* @brief MNMixtureClustering
*
* Fit a mixture of multivariate gaussian
* @author S. Riggi
* @date 30/07/2013
*/

#include <MNMixtureClustering.h>
#include <DataReader.h>
#include <MeanImputation.h>
#include <ListwiseDeletion.h>
#include <MultipleImputation.h>
#include <KMeansClustering.h>
#include <Util.h>
#include <MathUtils.h>

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

ClassImp(MDImputation_ns::MNMixtureClustering)

namespace MDImputation_ns {

//Static vars
//int MNMixtureClustering::fNComponents;
int MNMixtureClustering::fNDim;
long int MNMixtureClustering::fN;
//std::vector<TMatrixD> MNMixtureClustering::fData;
std::string MNMixtureClustering::fFileName;
MNClusteringOptions MNMixtureClustering::fOptions;

std::vector<double> MNMixtureClustering::fP;
std::vector<TMatrixD> MNMixtureClustering::fMu;
std::vector<TMatrixD> MNMixtureClustering::fSigma;
std::vector<TMatrixD> MNMixtureClustering::fSigmaInv;
std::vector<double> MNMixtureClustering::fSigmaDet;


MNMixtureClustering::MNMixtureClustering()
{
	fDataMatrix= 0;
	fDataMatrix_preImp= 0;
	fDataMatrix_imp= 0;
	fRTableName= "dataMatrix";
	fRTableName_preImp= "dataMatrix_preImp";
	fMeanData= 0;
	fCovarianceData= 0;
	fCovarianceDiagData= 0;
}

MNMixtureClustering::~MNMixtureClustering()
{
	//Clear data
	ClearData();
	
}//close destructor


void MNMixtureClustering::ClearData()
{
	//Clear data
	if(fDataMatrix){
		delete fDataMatrix;
		fDataMatrix= 0;
	}
	if(fDataMatrix_preImp){
		delete fDataMatrix_preImp;
		fDataMatrix_preImp= 0;
	}
	if(fDataMatrix_imp){
		delete fDataMatrix_imp;
		fDataMatrix_imp= 0;
	}
	if(fMeanData){
		delete fMeanData;
		fMeanData= 0;
	}
	if(fCovarianceData){
		delete fCovarianceData;
		fCovarianceData= 0;
	}
	if(fCovarianceDiagData){
		delete fCovarianceDiagData;
		fCovarianceDiagData= 0;
	}

}//close ClearData()


int MNMixtureClustering::Init()
{
	//## Clear R environment
	cout<<"INFO: Clearing R environment..."<<endl;
	Util::ClearRData();

	//### Initialize all necessary libraries
	cout<<"INFO: Loading needed R libraries..."<<endl;
	//std::vector<std::string> RLibraries{"mvtnorm","tmvtnorm","fMultivar","Matrix","moments"};
	std::vector<std::string> RLibraries{"Matrix"};
	if(Util::LoadRLibraries(RLibraries)<0){
		cerr<<"ERROR: Failed to load one or more of these R libraries {mvtnorm,tmvtnorm,fMultivar,Matrix,moments}, check if they are installed!"<<endl;
		return -1;
	}

	//### Initialize random generator
	delete gRandom;
	gRandom= new TRandom3(0);

	return 0;

}//close Init()

int MNMixtureClustering::ReadData()
{
	//## Read ascii data to ROOT matrix
	cout<<"INFO: Start reading data ..."<<endl;
	fDataMatrix= DataReader::ReadAscii(fFileName,fOptions.dataDelimiter);
	if(!fDataMatrix){
		cerr<<"ERROR: Failed to read data from file "<<fFileName<<"!"<<endl;
		return -1;
	}
	fNDim= fDataMatrix->GetNcols();
	fN= fDataMatrix->GetNrows();
	cout<<"INFO: Read "<<fN<<" x "<<fNDim<<" data table..."<<endl;

	//## Create matrix data copy (to be filled with imputed values)
	fDataMatrix_imp= new TMatrixD(fN,fNDim);

	//## Store index with observed and missing patterns in data matrix	
	cout<<"INFO: Store indexes with observed and missing patterns..."<<endl;
	for(long int i=0;i<fN;i++){	
		fObsDataIndexList.push_back( std::vector<long int>() );
		fMissingDataIndexList.push_back( std::vector<long int>() );

		//TMatrixD dataMatrix(fNDim,1);
		TMatrixD dataMatrix(1,fNDim);
		dataMatrix.Zero();

		for(int j=0;j<fNDim;j++){
			double dataValue= (*fDataMatrix)(i,j);
			(*fDataMatrix_imp)(i,j)= dataValue;
			//dataMatrix(j,0)= dataValue;
			dataMatrix(0,j)= dataValue;
			if(TMath::IsNaN(dataValue)){
				fMissingDataIndexList[i].push_back(j);	
			}
			else{
				fObsDataIndexList[i].push_back(j);	
			}
		}//end loop dim
		
		fData.push_back(dataMatrix);
		fData_completed.push_back(dataMatrix);
		
		//Store matrix with missing & obs info
		int nDim_obs= (int)(fObsDataIndexList[i].size());
		int nDim_miss= (int)(fMissingDataIndexList[i].size());

		//TMatrixD obsDataMatrix(nDim_obs,1);
		TMatrixD obsDataMatrix(1,nDim_obs);
		TMatrixD Sigma_mm(nDim_miss,nDim_miss);
		Sigma_mm.Zero();

		for(int s=0;s<nDim_obs;s++){	
			int VariableId_obs= fObsDataIndexList[i][s];
			//obsDataMatrix(s,0)= fData[i](VariableId_obs,0);
			obsDataMatrix(0,s)= fData[i](0,VariableId_obs);
		}//end loop observed dim

		fData_obs.push_back(obsDataMatrix);
		fProdSigmaMiss.push_back( std::vector<TMatrixD>() );
		fData_compl.push_back( std::vector<TMatrixD>() );

		for(int k=0;k<fOptions.nComponents;k++){
			fProdSigmaMiss[i].push_back(Sigma_mm);
			fData_compl[i].push_back(dataMatrix);			
		}
		
	}//end loop events

	//## Import matrix also in R for later use
	cout<<"INFO: Reading data table in R..."<<endl;
	if(DataReader::ReadAsciiInR(fFileName,fOptions.dataDelimiter,fRTableName)<0){
		cerr<<"ERROR: Failed to read ascii file and import it as an R table with name "<<fRTableName<<"!"<<endl;
		return -1;
	}

	//## Complete missing data with pre-imputation method
	cout<<"INFO: Running pre-imputation method "<<fOptions.preImputationMethod<<" ..."<<endl;
	fDataMatrix_preImp= nullptr;
	if(fOptions.preImputationMethod==MNClusteringOptions::eMEAN) {
		fDataMatrix_preImp= MeanImputation::RunImputationFromRTable(fRTableName,fRTableName_preImp);
	}
	else if(fOptions.preImputationMethod==MNClusteringOptions::eLD)	{
		fDataMatrix_preImp= ListwiseDeletion::RunImputationFromRTable(fRTableName,fRTableName_preImp);
	}
	else {
		cerr<<"ERROR: Invalid pre-impute method "<<fOptions.preImputationMethod<<"!"<<endl;
		return -1;
	}
	if(!fDataMatrix_preImp){
		cerr<<"ERROR: Failed to pre-impute data with method "<<fOptions.preImputationMethod<<"!"<<endl;
		return -1;
	}

	//## Compute covariance matrix & means of completed data with pre-imputation method
	fMeanData= Util::ComputeRTableColMeans(fRTableName_preImp);
	if(!fMeanData){
		cerr<<"ERROR: Failed to compute means of pre-imputed data!"<<endl;
		return -1;
	}
	cout<<"*** PRE-IMPUTED DATA - MEANS ***"<<endl;
	fMeanData->Print();
	
	fCovarianceData= Util::ComputeCovarianceMatrixFromRTable(fRTableName_preImp);
	if(!fCovarianceData){
		cerr<<"ERROR: Failed to compute covariance matrix of pre-imputed data!"<<endl;
		return -1;
	}
	cout<<"*** PRE-IMPUTED DATA - COV MATRIX ***"<<endl;
	fCovarianceData->Print();

	return 0;

}//close ReadData()


int MNMixtureClustering::RunImputation(std::string filename,MNClusteringOptions options)
{
	//Set options
	fFileName= filename;
	fOptions= options;

	//## Init data structures
	if(Init()<0){
		cerr<<"ERROR: Initialization failed!"<<endl;
		return -1;
	}
	
	//## Read input data
	if(ReadData()<0){
		cerr<<"ERROR: Failed to read data from file "<<fFileName<<"!"<<endl;
		return -1;
	}

	//## Run the EM algorithm
	if(RunEMClustering()<0){
		cerr<<"ERROR: EM clustering stage failed!"<<endl;
		return -1;
	}
	
	/*
	//## Compute classification info
	Classification();

	//## Draw fit info
	Draw();

	//## Save fit info
	Save();
	*/

	return 0;

}//close RunImputation()




int MNMixtureClustering::RunEMClustering()
{
	//## EM initialization phase
	if(RunEM_Init()<0){
		cerr<<"ERROR: EM initialization step failed!"<<endl;
		return -1;
	}

	/*
	//## Compute matrix with observed and missing data
	if(fUseMissingDataEM) ComputeObsMissData();
	
	//## Start EM iteration loop
	Print();
	*/

	//## Start EM iteration loop
	cout<<"INFO: Starting EM iteration loop (#"<<fOptions.nIterations<<" max niters) ..."<<endl;
	fLogLikelihood= 0.;

	for(int iter=0;iter<fOptions.nIterations;iter++)
	{
		//if(fUseTruncatedEM) ComputeGaussianTruncCorrection();

		//############################################################
		//## E step: Compute the expectation e0=tau, LL
		//############################################################			
		double LL= 0;
		if(RunEM_EStep(LL)<0){
			cerr<<"ERROR: EM EStep failed at iter no. "<<iter+1<<"!"<<endl;
			return -1;
		}
		double DeltaLogL= LL-fLogLikelihood;
	
		if(iter>1 && fOptions.useStoppingCriteria && DeltaLogL<fOptions.epsilon){
			cout<<"INFO: Stop criteria matched (LL="<<LL<<" DeltaLogL="<<DeltaLogL<<"< eps="<<fOptions.epsilon<<")...exit iteration!"<<endl;
			break;
		}

		fLogLikelihood= LL;
		fIterLogLikelihood.push_back(fLogLikelihood);
		//fIterNo= iter;
		//fFitIterInfo->Fill();
		//fIterLogLikelihood[iter]= fLogLikelihood;

		//############################################################
		//## M step: Update the fit parameters
		//############################################################
		RunEM_MStep();
		
		//############################################################
		//## Constrain Step
		//############################################################
		if(fOptions.useConstraints && RunEM_ConstrainStep()<0){
			cerr<<"ERROR: Failed to run EM constraint step!"<<endl;
			return -1;
		}

		//############################################################
		//## DUMP INFO
		//############################################################			
		cout<<"*** ITER NO. "<<iter<<" ***"<<endl;		
		PrintPars2();
		cout<<"==> LL= "<<LL<<"  DeltaL="<<DeltaLogL<<endl;
		cout<<"***************************"<<endl;
		cout<<endl;

	}//end loop EM iterations

	//## Correct the component weights with truncated EM
	/*
	if(fUseTruncatedEM){

		cout<<"p(";
		double normFactor= 0;
		for(int k=0;k<fNComponents;k++){
			TMatrixD a(fNDim,1);
			a= *fDataCutMin;
			TMatrixD b(fNDim,1);
			b= *fDataCutMax;
			
			//Compute the cumulative gaussian integrals
			double phi= MathUtilities::GetTruncatedGaussianCDF(fMu[k],fSigma[k],a,b);		
			double eta= fP[k]; 
			cout<<fP[k]<<",";
			normFactor+= eta/phi;
			fP[k]= eta/phi;
		}//end loop components
		cout<<")"<<endl;

		cout<<"p corr(";
		for(int k=0;k<fNComponents;k++) {
			fP[k]/= normFactor;
			cout<<fP[k]<<",";	
		}
		cout<<")"<<endl;

	}//close if
	*/

	
	//## Fill missing values in data with missing value estimates
	for(long int i=0;i<fN;i++)
	{
		int nDim_obs= (int)(fObsDataIndexList[i].size());
		int nDim_miss= (int)(fMissingDataIndexList[i].size());
				
		if(nDim_obs<=0){
			cout<<"--> No observed data...skip event!"<<endl;
			continue;
		}

		if(nDim_miss<=0) continue;
	
		for(int s=0;s<nDim_miss;s++){
			double missValueRec= 0.;

			for(int k=0;k<fOptions.nComponents;k++){
				double tau= fTau[i][k];
				//double Xm= fData_compl[i][k](fMissingDataIndexList[i][s],0);
				double Xm= fData_compl[i][k](0,fMissingDataIndexList[i][s]);
				//double weightedXm= fP[k]*Xm;//WRONG
				double weightedXm= tau*Xm;
				missValueRec+= weightedXm;
			}//end loop components
			//fData_completed[i](fMissingDataIndexList[i][s],0)= missValueRec;
			fData_completed[i](0,fMissingDataIndexList[i][s])= missValueRec;
				
		}//end loop miss dim
	
		//for(int j=0;j<fNDim;j++) (*fDataMatrix_imp)(i,j)= fData_completed[i](j,0);
		for(int j=0;j<fNDim;j++) (*fDataMatrix_imp)(i,j)= fData_completed[i](0,j);
	
	}//end loop events
	
	return 0;

}//close RunEMClustering()


int MNMixtureClustering::RunEM_Init()
{	
	//## Init vars
	for(long int i=0;i<fN;i++){
		fTau.push_back( std::vector<double>() );
		
		for(int k=0;k<fOptions.nComponents;k++){
			fTau[i].push_back(0.);
		}//end loop components
	}//end loop data  
	
	
	//## Init start mixture parameters
	for(int k=0;k<fOptions.nComponents;k++){
		//TMatrixD startMu(fNDim,1);
		TMatrixD startMu(1,fNDim);
		startMu.Zero();

		TMatrixD startSigma(fNDim,fNDim);		
		startSigma.Zero();

		//Initialize truncated gaussian correction pars
		fTruncGaussianNormFactor.push_back(0.);
		fMuTruncGaussianCorrection.push_back(startMu);
		fSigmaTruncGaussianCorrection.push_back(startSigma);
		//fInitRandNumbers.push_back(0);

		//Initialize starting pars
		fMu_start.push_back(startMu);
		fSigma_start.push_back(startSigma);
		fSigmaInv_start.push_back(startSigma);
		fP_start.push_back(1./fOptions.nComponents);

		//Initialize starting eigen
		fSigmaEigen_start.push_back(startMu);
		fSigmaEigenvect_start.push_back(startSigma);

		/*
		fMu_true.push_back(startMu);
		fSigma_true.push_back(startSigma);
		fP_true.push_back(1./fNComponents);
		*/

		//Initialize fitted pars
		fMu.push_back(startMu);
		fMuDiff.push_back(startMu);	
		fSigma.push_back(startSigma);
		fSigmaInv.push_back(startSigma);
		fSigmaDet.push_back(0.);
		fSigmaDiff.push_back(startSigma);	
		fSigmaEigen.push_back(startMu);
		fSigmaEigenvect.push_back(startSigma);
		fP.push_back(1./fOptions.nComponents);
		
		//Initialize safe par values
		fMu_safe.push_back(startMu);
		fSigma_safe.push_back(startSigma);	
		fP_safe.push_back(1./fOptions.nComponents);
		fSigmaEigen_safe.push_back(startMu);
	
	}//end loop components

	/*
	//## Set start parameters
	fP= fP_startFile;
	fMu= fMu_startFile;
	fSigma= fSigma_startFile;
	*/

	//## Initialize mixture pars
	int status= 0;
	if(fOptions.parInitMethod==MNClusteringOptions::eKMEANS){
		status= InitParsToKMeans();
	}
	else if(fOptions.parInitMethod==MNClusteringOptions::eRANDOM){	
		cerr<<"ERROR: Randomized parameter initialization not yet implemented!"<<endl;
		status= -1;		
	}
	else if(fOptions.parInitMethod==MNClusteringOptions::eUSER){
		//cerr<<"ERROR: User-provided parameter initialization not yet implemented!"<<endl;
		status= InitParsToUser();
	}
	else{
		cerr<<"ERROR: Invalid par initialization method given ("<<fOptions.parInitMethod<<")!"<<endl;
		return -1;
	}

	if(status<0){
		cerr<<"ERROR: EM parameter initialization failed!"<<endl;
		return -1;
	}

	//## Randomize pars?
	if(fOptions.randomizeStartPars){
		RandomizePars();
	}	

	/*
	//## Use kmeans initialization?
	if(fUseKMeansStartPar){
		InitParsToKMeans();
	}
	
	//## Randomize start parameters
	if(fRandomizeStartPar){
		RandomizePar();
	}

	//## Randomize start parameters from models
	if(fUseRandomFromModelStart){
		RandomizeParFromModel();
	}
	*/
	
	//## Compute sigma inverse & eigenvalues
	for(int k=0;k<fOptions.nComponents;k++) 
	{	
		//Compute inverse
		fSigmaDet[k]= fSigma[k].Determinant();	
		fSigmaInv[k]= TMatrixD(TMatrixD::kInverted,fSigma[k]);
		if (fSigmaDet[k]<=0) {
			cerr<<"WARN: Covariance matrix inversion failed for component "<<k+1<<" (SigmaDet="<<fSigmaDet[k]<<")!"<<endl;
			return -1;
		}	

		//Compute eigenvalues
		//NB: Fails if imput matrix not symmetric
		if(Util::ComputeSymMatrixEigenvalues(fSigmaEigen[k],fSigmaEigenvect[k],fSigma[k])<0){
			cerr<<"ERROR: Failed to compute covariance matrix eigenvalues for component "<<k+1<<"!"<<endl;
			return -1;
		}

		//## Check if sigma eigen fullfil the constraints
		if(fOptions.useConstraints && fOptions.useCovarianceEigenBoundConstraint){
			for(int j=0;j<fNDim;j++){
				double lambda= fSigmaEigen[k](0,j);
				double lambdaMin= fOptions.SigmaEigen_min[k](0,j);
				double lambdaMax= fOptions.SigmaEigen_max[k](0,j);

				if(lambda<=lambdaMin || lambda>=lambdaMax){
					cerr<<"ERROR: Eigenvalues for component no. "<<k+1<<" not satisfying the constraints for component "<<k+1<<" (lambda="<<lambda<<" min/max="<<lambdaMin<<"/"<<lambdaMax<<")...exit!"<<endl;
					return -1;
				}
			}//end loop ndim
		}//close if
		
	}//end loop components

	/*
	//## Check if mean and covariance matrix fullfil the constraints
	if(fUseCovarianceConstraint){

		for(int j=0;j<fNDim;j++){
			double constraintSign= 1;
			
			for(int k=0;k<fNComponents-1;k++){
				double Sigma_k= fSigma[k](j,j);
				double Sigma_k_1= fSigma[k+1](j,j);
				
				if(constraintSign*Sigma_k<=constraintSign*Sigma_k_1) {
					cerr<<"MGMixtureFitter::EMInit(): ERROR: Covariance constraint not valid at start for component "<<k+1<<" (Sigma(k)="<<Sigma_k<<" Sigma(k+1)="<<Sigma_k_1<<") ...exit!"<<endl;
					exit(1);	
				}
					
			}//end loop components	
		}//end loop dim

	}//close if	

	//## Check if covariance fullfil the constraints
	if(fUseCovarianceBoundConstraint){

		for(int k=0;k<fNComponents;k++){

			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++){
					double constraintSign= 1;
			
					double Sigma= fSigma[k](j,l);
					double Sigma_min= fSigma_min[k](j,l);
					double Sigma_max= fSigma_max[k](j,l);
				
					
					if(constraintSign*Sigma<=constraintSign*Sigma_min) {
						cerr<<"MGMixtureFitter::EMInit(): ERROR: Covariance min bound constraint not valid at start for component "<<k+1<<" (Sigma="<<Sigma<<" Sigma_min="<<Sigma_min<<") ...exit!"<<endl;
						exit(1);	
					}
					if(constraintSign*Sigma>=constraintSign*Sigma_max) {
						cerr<<"MGMixtureFitter::EMInit(): ERROR: Covariance max bound constraint not valid at start for component "<<k+1<<" (Sigma="<<Sigma<<" Sigma_max="<<Sigma_max<<") ...exit!"<<endl;
						exit(1);	
					}
					
				}//end loop dim
			}//end loop dim
		}//end loop components	
	}//close if

	//## Check if mean vector fullfil the constraints
	if(fUseMeanConstraint){

		for(int j=0;j<fNDim;j++){
			
			double constraintSign= (*fMeanConstraintSign)(j,0);
 
			for(int k=0;k<fNComponents-1;k++){
				double Mu_k= fMu[k](j,0);
				double Mu_k_1= fMu[k+1](j,0);
				
				if(constraintSign*Mu_k<=constraintSign*Mu_k_1) {
					cerr<<"MGMixtureFitter::EMInit(): ERROR: Mean constraint not valid at start (Mu(k)="<<Mu_k<<" Mu(k+1)="<<Mu_k_1<<") ...exit!"<<endl;
					exit(1);	
				}
					
			}//end loop components	
		}//end loop dim
	}//close if

	//## Check if mean vector fullfil the constraints
	if(fUseMeanBoundConstraint){

		for(int k=0;k<fNComponents;k++){

			for(int j=0;j<fNDim;j++){
				double constraintSign= 1;
			
				double Mu= fMu[k](j,0);
				double Mu_min= fMu_min[k](j,0);
				double Mu_max= fMu_max[k](j,0);
				
				if(constraintSign*Mu<=constraintSign*Mu_min) {
					cerr<<"MGMixtureFitter::EMInit(): ERROR: Mean min bound constraint not valid at start (Mu="<<Mu<<" Mu_min="<<Mu_min<<") ...exit!"<<endl;
					exit(1);	
				}
				if(constraintSign*Mu>=constraintSign*Mu_max) {
					cerr<<"MGMixtureFitter::EMInit(): ERROR: Mean max bound constraint not valid at start (Mu="<<Mu<<" Mu_max="<<Mu_max<<") ...exit!"<<endl;
					exit(1);	
				}
					
			} //end loop dim
		}//end loop components	
	}//close if
	*/

	//## Assign starting values
	cout<<"MSTMixtureFitter::EMInit(): Starting parameters..."<<endl;
	for(int k=0;k<fOptions.nComponents;k++){	
		
		fMu_start[k]= fMu[k];
		fSigma_start[k]= fSigma[k];
		fP_start[k]= fP[k];
		fSigmaInv_start[k]= fSigmaInv[k];
		fSigmaEigen_start[k]= fSigmaEigen[k];
		fSigmaEigenvect_start[k]= fSigmaEigenvect[k];

		cout<<"== Component "<<k+1<<" =="<<endl;
		cout<<"p= "<<fP[k]<<endl;
		cout<<"Mu= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(0,j)<<",";
		cout<<(fMu[k])(0,fNDim-1)<<")"<<endl; 

		cout<<"Sigma= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
		cout<<")"<<endl;

		cout<<"SigmaInv= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigmaInv[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
		
		cout<<"SigmaEigen= (";
		for(int j=0;j<fNDim;j++){
			cout<<(fSigmaEigen[k])(0,j)<<",";
		}	
		cout<<")"<<endl;
	}

	//## Set safe values
	for(int k=0;k<fOptions.nComponents;k++){
		fMu_safe[k]= fMu[k];
		fSigma_safe[k]= fSigma[k];
		fSigmaEigen_safe[k]= fSigmaEigen[k];	
		fP_safe[k]= fP[k];
	}	

	return 0;

}//close RunEM_Init()


int MNMixtureClustering::RunEM_EStep(double& LL)
{
	//## Compute inverse of covariance matrix
	for(int k=0;k<fOptions.nComponents;k++){
		double SigmaDet= 0.;
		
		TMatrixD SigmaInv(fNDim,fNDim);
		TMatrixD Sigma(fNDim,fNDim);
		Sigma= fSigma[k];

		SigmaInv = Sigma.Invert(&SigmaDet);
			
		fSigmaInv[k]= SigmaInv;		
		fSigmaDet[k]= SigmaDet;	
  		
		if (SigmaDet<=0) {
			cerr<<"WARN: Covariance matrix inversion failed for component "<<k+1<<" Sigma=(";
			for(int l=0;l<fNDim;l++) {
				for(int j=0;j<fNDim;j++) cout<<fSigma[k](l,j)<<",";
			}
			cerr<<")  SigmaDet="<<fSigmaDet[k]<<")!"<<endl;
			return -1;
		}		
	}//end loop mixtures


	//## Compute posterior probability & LL
	LL= 0.;
	for(long int i=0;i<fN;i++){
			
		double tauSum= 0.;

		//## Check for existing observed data
		int nDim_obs= (int)(fObsDataIndexList[i].size());
		int nDim_miss= (int)(fMissingDataIndexList[i].size());
		if(nDim_obs<=0) {
			cout<<"--> No observed data...skip event!"<<endl;
			continue;
		}
				
		//## Loop over components
		for(int k=0;k<fOptions.nComponents;k++){

			double tauComponent= 0.;

			//## In case of missing data compute for each event the observed mu, sigma, sigmaInv matrix and p
			//TMatrixD Mu_obs(nDim_obs,1);
			TMatrixD Mu_obs(1,nDim_obs);
			Mu_obs.Zero();
			TMatrixD Sigma_obs(nDim_obs,nDim_obs);
			Sigma_obs.Zero();
			TMatrixD SigmaInv_obs(nDim_obs,nDim_obs);
			SigmaInv_obs.Zero();
			double SigmaDet_obs= 0.;
	
			if(nDim_miss>0){
				for(int s=0;s<nDim_obs;s++){
					//Mu_obs(s,0)= fMu[k](fObsDataIndexList[i][s],0); 	
					Mu_obs(0,s)= fMu[k](0,fObsDataIndexList[i][s]); 	
					
					for(int t=0;t<nDim_obs;t++){
						Sigma_obs(s,t)= fSigma[k](fObsDataIndexList[i][s],fObsDataIndexList[i][t]);
					}//end loop obs dim
				}//end loop obs dim

				SigmaInv_obs= TMatrixD(TMatrixD::kInverted, Sigma_obs);
				SigmaDet_obs= Sigma_obs.Determinant();
				//SigmaInv_obs= Sigma_obs;
				//SigmaInv_obs= SigmaInv_obs.Invert(&SigmaDet_obs);
				if (SigmaDet_obs<=0) {
					cerr<<"WARN: Observed Covariance matrix inversion failed!"<<endl;
				}

				tauComponent= MathUtils::tauGaus(fData_obs[i],Mu_obs,SigmaInv_obs,SigmaDet_obs,fP[k]);
				
			}//close if have missing data 
			else{
				tauComponent= MathUtils::tauGaus(fData[i],fMu[k],fSigmaInv[k],fSigmaDet[k],fP[k]);
			}
		
			//if(fUseTruncatedEM){//check if truncated EM is used
			//	fTau[i][k]= tauComponent/fTruncGaussianNormFactor[k];
			//	tauSum+= tauComponent/fTruncGaussianNormFactor[k];	
			//}
			//else{//standard EM computation
				fTau[i][k]= tauComponent;
				tauSum+= tauComponent;	
			//}
						
		}//end loop components

		//## Compute LL
		LL+= log(tauSum);
			
		for(int k=0;k<fOptions.nComponents;k++) {
			if(tauSum>0) fTau[i][k]/= tauSum;
			else{
				cerr<<"ERROR: Negative or null tau sum for event no. "<<i<<" (tauSum="<<tauSum<<" data(";
				for(int j=0;j<fNDim;j++){
					//cerr<<fData[i](j,0)<<",";	
					cerr<<fData[i](0,j)<<",";
				}
				cerr<<")...exit!"<<endl;
				return -1;
			}
		}//end loop mixture components
		
	}//end loop events

	return 0;

}//close RunEM_EStep()


int MNMixtureClustering::RunEM_MStep()
{
	//## Compute the EM parameter update
	for(int k=0;k<fOptions.nComponents;k++)
	{	
		double tauSum= 0.;
		//TMatrixD muSum(fNDim,1);
		TMatrixD muSum(1,fNDim);
		muSum.Zero();
		TMatrixD sigmaSum(fNDim,fNDim);
		sigmaSum.Zero();

		TMatrixD diff(1,fNDim);
		diff.Zero();
		TMatrixD diff_t(fNDim,1);
		diff_t.Zero();
		
		/*
		TMatrixD diff(fNDim,1);
		diff.Zero();
		TMatrixD diff_t(1,fNDim);
		diff_t.Zero();
		*/

		//## Loop over events	
		for(long int i=0;i<fN;i++){
			int nDim_obs= (int)(fObsDataIndexList[i].size());
			int nDim_miss= (int)(fMissingDataIndexList[i].size());
				
			if(nDim_obs<=0){
				cout<<"--> No observed data...skip event!"<<endl;
				continue;
			}

			//## Use missing data EM?
			if(nDim_miss>0){
				TMatrixD Data_obs(nDim_obs,1);
				TMatrixD Data_miss(nDim_miss,1);
				TMatrixD Mu_miss(nDim_miss,1);
				TMatrixD Mu_obs(nDim_obs,1);
				
				TMatrixD Sigma_mm(nDim_miss,nDim_miss);
				TMatrixD Sigma_mo(nDim_miss,nDim_obs);
				TMatrixD Sigma_mo_t(nDim_obs,nDim_miss);
				TMatrixD Sigma_oo(nDim_obs,nDim_obs);
				TMatrixD SigmaInv_oo(nDim_obs,nDim_obs);
				TMatrixD SigmaMissProd(nDim_miss,nDim_miss);
	
				double SigmaDet_oo= 0.;

				for(int s=0;s<nDim_obs;s++){
					Data_obs(s,0)= fData_obs[i](0,s);

					//Mu_obs(s,0)= fMu[k](fObsDataIndexList[i][s],0); 
					Mu_obs(s,0)= fMu[k](0,fObsDataIndexList[i][s]); 

					for(int t=0;t<nDim_obs;t++){
						Sigma_oo(s,t)= fSigma[k](fObsDataIndexList[i][s],fObsDataIndexList[i][t]);
					}//end loop obs dim
	
					for(int t=0;t<nDim_miss;t++){
						Sigma_mo(t,s)= fSigma[k](fMissingDataIndexList[i][t],fObsDataIndexList[i][s]);
					}//end loop missing dim
				}//end loop obs dim

				for(int s=0;s<nDim_miss;s++){
					//Mu_miss(s,0)= fMu[k](fMissingDataIndexList[i][s],0); 
					Mu_miss(s,0)= fMu[k](0,fMissingDataIndexList[i][s]); 

					for(int t=0;t<nDim_miss;t++){
						Sigma_mm(s,t)= fSigma[k](fMissingDataIndexList[i][s],fMissingDataIndexList[i][t]);
					}//end loop missing dim
				}//end loop missing dim

				Sigma_mo_t= TMatrixD(TMatrixD::kTransposed, Sigma_mo);
				
				SigmaInv_oo= TMatrixD(TMatrixD::kInverted, Sigma_oo);
				SigmaDet_oo= Sigma_oo.Determinant();
				if(SigmaDet_oo<=0){
					cerr<<"WARN: Inversion of observed covariance failed!"<<endl;
				}	

				SigmaMissProd= Sigma_mm-Sigma_mo*SigmaInv_oo*Sigma_mo_t;
				fProdSigmaMiss[i][k]= fTau[i][k]*SigmaMissProd;

				//cout<<"INFO: Computing Data_miss (Data_obs="<<fData_obs[i].GetNrows()<<" x "<<fData_obs[i].GetNcols()<<", Mu_obs="<<Mu_obs.GetNrows()<<" x "<<Mu_obs.GetNcols()<<", SigmaInv_oo="<<SigmaInv_oo.GetNrows()<<" x "<<SigmaInv_oo.GetNcols()<<", Sigma_mo="<<Sigma_mo.GetNrows()<<" x "<<Sigma_mo.GetNcols()<<", Mu_miss="<<Mu_miss.GetNrows()<<" x "<<Mu_miss.GetNcols()<<")..."<<endl;
				//Data_miss= Mu_miss + Sigma_mo*SigmaInv_oo*(fData_obs[i]-Mu_obs);
				Data_miss= Mu_miss + Sigma_mo*SigmaInv_oo*(Data_obs-Mu_obs);
				//cout<<"Data_miss size= "<<Data_miss.GetNrows()<<" x "<<Data_miss.GetNcols()<<endl;

				for(int s=0;s<nDim_miss;s++){
					//fData_compl[i][k](fMissingDataIndexList[i][s],0)= Data_miss(s,0);
					fData_compl[i][k](0,fMissingDataIndexList[i][s])= Data_miss(s,0);
				}//end loop missing dim

				tauSum+= fTau[i][k];
				muSum+= fTau[i][k]*fData_compl[i][k];

			}//close if missing data
			else{//standard EM
				tauSum+= fTau[i][k];
				muSum+= fTau[i][k]*fData[i];
			}
		
		}//end loop events

		double tauSumInv= 1./tauSum;		
	
		//## Update fraction parameters
		if(!fOptions.fixFractionPars) fP[k]= tauSum/(double)(fN);

		//## Update mean parameters
		if(!fOptions.fixMeanPars) {
			fMu[k]= muSum*tauSumInv;

			//## Apply correction term for truncated EM
			//if(fUseTruncatedEM) fMu[k]= fMu[k] - fMuTruncGaussianCorrection[k];
					
		}//close !fFixMeanPar

		
		//## Update sigma 	
		if(!fOptions.fixCovariancePars){
			cout<<"INFO: Updating covariance pars..."<<endl;
			for(int i=0;i<fN;i++){
				diff= (fData_compl[i][k]-fMu[k]);			
				diff_t= TMatrixD(TMatrixD::kTransposed, diff);
				//sigmaSum+= fTau[i][k]*diff*diff_t;
				sigmaSum+= fTau[i][k]*diff_t*diff;
			}//end loop events
			
			fSigma[k]= sigmaSum*tauSumInv;
				
			//## Add covariance of missing data to the missing part of sigma	
			cout<<"INFO: Add covariance of missing data to the missing part of sigma	..."<<endl;
			for(long int i=0;i<fN;i++){
				int nDim_miss= (int)(fMissingDataIndexList[i].size());
				for(int s=0;s<nDim_miss;s++){
					for(int t=0;t<nDim_miss;t++){
						fSigma[k](fMissingDataIndexList[i][s],fMissingDataIndexList[i][t])+= fProdSigmaMiss[i][k](s,t)/tauSum;
					}//end loop missing dim
				}//end loop missing dim		
			}//end loop events
	
		}//close if !fFixCovariancePar
	}//end loop components

	//## Check covariance matrix integrity
	if(!fOptions.fixCovariancePars && CheckCovariance()<0){
		cerr<<"ERROR: Failed to perform checks in covariance matrix!"<<endl;
		return -1;
	}

	return 0;

}//close RunEM_MStep()


int MNMixtureClustering::RunEM_ConstrainStep()
{
	//==================================================
	//==      SIGMA EIGENVALUES BOUND CONSTRAINT
	//==================================================
	//## Check covariance eigenvalues constraints
	//## Eigenvalues must be in desired bounds
	if(fOptions.useCovarianceEigenBoundConstraint) {
		for(int k=0;k<fOptions.nComponents;k++)
		{
			bool isSigmaEigenBoundConstraintViolated= false;
			double alphaLambda_minmax[fNDim];
	
			for(int j=0;j<fNDim;j++){
				alphaLambda_minmax[j]= 1;

				double Lambda= fSigmaEigen[k](0,j);
				double Lambda_safe= fSigmaEigen_safe[k](0,j);
				double Lambda_min= fOptions.SigmaEigen_min[k](0,j);
				double Lambda_max= fOptions.SigmaEigen_max[k](0,j);
				double alphaLambda_minBound= (Lambda_min-Lambda_safe)/(Lambda-Lambda_safe);
				double alphaLambda_maxBound= (Lambda_max-Lambda_safe)/(Lambda-Lambda_safe);
			
				if(Lambda<Lambda_min){//min constraint violated
					alphaLambda_minmax[j]= alphaLambda_minBound;
					isSigmaEigenBoundConstraintViolated= true;
					cout<<"DEBUG: Covariance min eigen bound constraint violated for component no. "<<k+1<<":  eigen="<<Lambda<<"  eigenMin="<<Lambda_min<<"  eigenSafe="<<Lambda_safe<<"  alpha="<<alphaLambda_minBound<<endl;	
				}
				if(Lambda>Lambda_max){//max constraint violated
					alphaLambda_minmax[j]= alphaLambda_maxBound;
					isSigmaEigenBoundConstraintViolated= true;
					cout<<"DEBUG: Covariance max eigen bound constraint violated for component no. "<<k+1<<":  eigen="<<Lambda<<"  eigenMax="<<Lambda_max<<"  eigenSafe="<<Lambda_safe<<"  alpha="<<alphaLambda_maxBound<<endl;		
				}
			}//end loop dim

			//## Update covariance eigen values		
			if(isSigmaEigenBoundConstraintViolated){
			
				//## Find mix alpha value and define optimal alpha
				double alphaEigenMin= 1;
				cout<<"DEBUG: alphaSigmaEigenList (";
				for(int j=0;j<fNDim;j++){
					cout<<alphaLambda_minmax[j]<<",";
					if(alphaLambda_minmax[j]<alphaEigenMin) alphaEigenMin= alphaLambda_minmax[j];
				}				
				cout<<")"<<endl;
			
				double alphaEigenOpt= alphaEigenMin/fOptions.constraintAlphaScale;

				cout<<"DEBUG: alphaEigenMin="<<alphaEigenMin<<", alphaEigenOpt="<<alphaEigenOpt<<endl;

				//## Update sigma eigen values
				if(alphaEigenOpt<fOptions.constraintAlphaTolerance) {
					cout<<"WARN: Covariance eigen for component "<<k+1<<" is stuck in constraint ("<<alphaEigenOpt<<"<"<<fOptions.constraintAlphaTolerance<<")!"<<endl;		
					fSigmaEigen[k]= fSigmaEigen_start[k];
				}
				else{
					for(int j=0;j<fNDim;j++){
						double Lambda= fSigmaEigen[k](0,j);
						double Lambda_safe= fSigmaEigen_safe[k](0,j);
						fSigmaEigen[k](0,j)= (1.-alphaEigenOpt)*Lambda_safe + alphaEigenOpt*Lambda;
					}//end loop dim
				}
				
				//## Re-Calculate covariance starting from the updated eigenvalues
				fSigma[k]= fSigmaEigenvect[k]*fSigmaEigen[k]*TMatrixD(TMatrixD::kInverted,fSigmaEigenvect[k]);

			}//close if isSigmaEigenBoundConstraintViolated
		}//end loop components

		//## Check covariance (re-compute inverse, determinant, etc)
		if(CheckCovariance()<0){
			cerr<<"ERROR: Failed to check covariance after eigen constraint!"<<endl;
			return -1;
		}
				
	}//close fUseCovarianceEigenConstraint


	/*
	//###########################
	//####  MEAN CONSTRAINT #####
	//###########################
	//## Apply constraints on component means
	if(!fFixMeanPar){
	
		//### MEAN BOUND CONSTRAINT
		std::vector<TMatrixD> alphaMuBoundList;
		alphaMuBoundList.clear();
		alphaMuBoundList.resize(0);
		bool isMuBoundConstraintViolated= false;
		bool isBadMuBound[fNComponents];

		if(fUseMeanBoundConstraint){
		
			for(int k=0;k<fNComponents;k++){
				isBadMuBound[k]= false;

				TMatrixD alphaMuBound(fNDim,1);
				alphaMuBound.Zero();
	
				for(int j=0;j<fNDim;j++){
					alphaMuBound(j,0)= 1;
					double constraintSign= 1;
					
					double a_k= fMu_min[k](j,0);
					double b_k= fMu_max[k](j,0);
					
					double Mu= fMu[k](j,0);
					double Mu_safe= fMu_safe[k](j,0);
				
					double alphaMu_minBound= (a_k-Mu_safe)/(Mu-Mu_safe);
					double alphaMu_maxBound= (b_k-Mu_safe)/(Mu-Mu_safe);

					if(constraintSign*Mu<constraintSign*a_k){//min constraint violated
						alphaMuBound(j,0)= alphaMu_minBound;
						isBadMuBound[k]= true;
						isMuBoundConstraintViolated= true;
						cout<<"MGMixtureFitter::EMConstrain(): INFO: Mean min bound constraint violated for component no. "<<k+1<<":  Mu="<<Mu<<"  Mu_min="<<a_k<<"  Mu_safe="<<Mu_safe<<"  alpha="<<alphaMu_minBound<<endl;

						if(alphaMu_minBound<0 || alphaMu_minBound>1){
							cerr<<"MGMixtureFitter::EMConstrain(): ERROR: Invalid alpha value for mu min bound constraint for component no. "<<k+1<<" (alpha="<<alphaMu_minBound<<") ...exit!"<<endl;
							exit(1);
						}		
					}
					
					if(constraintSign*Mu>b_k){//max constraint violated
						alphaMuBound(j,0)= alphaMu_maxBound;
						isBadMuBound[k]= true;
						isMuBoundConstraintViolated= true;
						cout<<"MGMixtureFitter::EMConstrain(): INFO: Mean max bound constraint violated for component no. "<<k+1<<":  Mu="<<Mu<<"  Mu_max="<<b_k<<"  Mu_safe(k)="<<Mu_safe<<"  alpha="<<alphaMu_maxBound<<endl;	
						if(alphaMu_maxBound<0 || alphaMu_maxBound>1){
							cerr<<"MGMixtureFitter::EMConstrain(): ERROR: Invalid alpha value for mu max bound constraint for component no. "<<k+1<<" (alpha="<<alphaMu_maxBound<<") ...exit!"<<endl;
							exit(1);
						}	
					}	
					
				}//end loop dim

				alphaMuBoundList.push_back(alphaMuBound);
				
			}//end loop components
		}//close if use mean bound constraint


	
		//### MEAN GROUP CONSTRAINT
		std::vector<double> alphaMuList;
		alphaMuList.clear();
		alphaMuList.resize(0);
		bool isMuGroupConstraintViolated= false;

		//## Check if mean order constraint is violated
		if(fUseMeanConstraint){

			for(int j=0;j<fNDim;j++){
				//double constraintSign= 1;
				//if(j==1) constraintSign= -1;
				double constraintSign= (*fMeanConstraintSign)(j,0);

				for(int k=0;k<fNComponents-1;k++){
					double Mu_k= fMu[k](j,0);
					double Mu_k_1= fMu[k+1](j,0);
					double safeMu_k= fMu_safe[k](j,0);
					double safeMu_k_1= fMu_safe[k+1](j,0);
					double denom= 1.-(Mu_k_1-Mu_k)/(safeMu_k_1-safeMu_k);
					double alpha= 1./denom;
		
					if(constraintSign*Mu_k>constraintSign*Mu_k_1) continue;//constraint satisfied...skip
		
					cout<<"MGMixtureFitter::EMConstrain(): INFO: Mean constraint violated for component no. "<<k+1<<":  mu(k)="<<Mu_k<<"  mu(k+1)="<<Mu_k_1<<"  mu_safe(k)="<<safeMu_k<<"  mu_safe(k+1)="<<safeMu_k_1<<"  alpha="<<alpha<<endl;

					if(alpha<0 || alpha>1){
						cerr<<"MGMixtureFitter::EMConstrain(): ERROR: Invalid alpha value for mu group constraint for component "<<k+1<<"-"<<k+2<<" (alpha="<<alpha<<") ...exit!"<<endl;
						exit(1);
					}	

					isMuGroupConstraintViolated= true;
					alphaMuList.push_back(alpha);

				}//end loop components	
			}//end loop dim

		}//close if use mean constraint


		//## Computing min alpha
		double minAlpha_groupConstraint= 1;
		if(isMuGroupConstraintViolated){
			cout<<"MGMixtureFitter::EMConstrain(): INFO: alphaMuList(";
			for(unsigned int i=0;i<alphaMuList.size();i++){
				cout<<alphaMuList[i]<<",";
				if(alphaMuList[i]<minAlpha_groupConstraint) minAlpha_groupConstraint= alphaMuList[i];
			}	
			cout<<")"<<endl;
		}
		
		double minAlpha_boundConstraint[fNComponents];
		for(int k=0;k<fNComponents;k++){
			minAlpha_boundConstraint[k]= 1;
			if(isBadMuBound[k]){
				minAlpha_boundConstraint[k]= alphaMuBoundList[k].Min();	
			}
		}//end loop components

		//## Updating the mean vector
		if(isMuGroupConstraintViolated || isMuBoundConstraintViolated){
			cout<<"MGMixtureFitter::EMConstrain(): INFO: Updating mean vector..."<<endl;
		
			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateMu= false;

			for(int k=0;k<fNComponents;k++){
				if(!isMuGroupConstraintViolated && !isBadMuBound[k]) continue;

				double alphaMuMin= min(minAlpha_groupConstraint,minAlpha_boundConstraint[k]);
				double alphaMuOpt= alphaMuMin/fConstraintAlphaScale;

				cout<<"MGMixtureFitter::EMConstrain(): INFO: Component no. "<<k+1<<": minAlpha_groupConstraint="<<minAlpha_groupConstraint<<"  minAlpha_boundConstraint="<<minAlpha_boundConstraint[k]<<"  alphaMuMin="<<alphaMuMin<<"  alphaMuOpt="<<alphaMuOpt<<endl;

				if(alphaMuOpt<fConstraintAlphaTolerance) {
					cerr<<"MGMixtureFitter::EMConstrain(): WARNING: Mean vector for component "<<k+1<<" is stuck in constraint ("<<alphaMuOpt<<"<"<<fConstraintAlphaTolerance<<")... ...regenerate!"<<endl;
					hasToRegenerateMu= true;
					hasToBeGenerated[k]= true;	
				}
			
				for(int j=0;j<fNDim;j++){
					double Mu= fMu[k](j,0);
					double safeMu= fMu_safe[k](j,0);
					fMu[k](j,0)= (1.-alphaMuOpt)*safeMu + alphaMuOpt*Mu;
				}//end loop dim		
			}//end loop components
		
			if(hasToRegenerateMu && fUseRandomRegenerationAfterStuck){
				cout<<"MGMixtureFitter::EMConstrain(): INFO: Regenerating mu pars..."<<endl;
				RandomizeMeanParFromModel(hasToBeGenerated);	
			}

		}//close if constraint violated

	}//close if fFixMeanPar


	//############################
	//##   SIGMA CONSTRAINTS    ##
	//############################
	cout<<"Sigma Constraint"<<endl;
	//## Apply constraints on component covariances
	if(!fFixCovariancePar){

		
		bool isBadSigmaBound[fNComponents];
		for(int k=0;k<fNComponents;k++) isBadSigmaBound[k]= false;
		bool isSigmaBoundConstraintViolated= false;
		std::vector<TMatrixD> alphaSigmaBoundList;
		alphaSigmaBoundList.clear();
		alphaSigmaBoundList.resize(0);

		if(fUseCovarianceBoundConstraint){

			for(int k=0;k<fNComponents;k++){
				
				isBadSigmaBound[k]= false;

				TMatrixD alphaSigmaBound(fNDim,fNDim);
				alphaSigmaBound.Zero();

				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double constraintSign= 1;
						alphaSigmaBound(j,l)= 1;


						//if(j!=l) continue;//Do not apply constraint on covariance terms

						
						double Sigma= fSigma[k](j,l);
						double Sigma_safe= fSigma_safe[k](j,l);
						double Sigma_min= fSigma_min[k](j,l);
						double Sigma_max= fSigma_max[k](j,l);
						double alphaSigma_minBound= (Sigma_min-Sigma_safe)/(Sigma-Sigma_safe);
						double alphaSigma_maxBound= (Sigma_max-Sigma_safe)/(Sigma-Sigma_safe);
		
						
						if(constraintSign*Sigma<constraintSign*Sigma_min){//min constraint violated
							alphaSigmaBound(j,l)= alphaSigma_minBound;
							isBadSigmaBound[k]= true;
							isSigmaBoundConstraintViolated= true;

							if(j==l) cout<<"INFO: Sigma min constraint violated for component no. "<<k+1<<":  Sigma(k)="<<Sigma<<"  Sigma_min(k)="<<Sigma_min<<"  Sigma_safe(k)="<<Sigma_safe<<"  alpha="<<alphaSigma_minBound<<endl;	
							else cout<<"INFO: Covariance min constraint violated for component no. "<<k+1<<":  Sigma(k)="<<Sigma<<"  Sigma_min(k)="<<Sigma_min<<"  Sigma_safe(k)="<<Sigma_safe<<"  alpha="<<alphaSigma_minBound<<endl;	

							if(alphaSigma_minBound<0 || alphaSigma_minBound>1){
								cerr<<"MSNMixtureFitter::EMConstrain(): ERROR: Invalid alpha value for covariance min bound constraint for component "<<k+1<<"-"<<k+2<<" (alpha="<<alphaSigma_minBound<<") ...exit!"<<endl;
								exit(1);
							}

						}
						if(constraintSign*Sigma>constraintSign*Sigma_max){//max constraint violated
							alphaSigmaBound(j,l)= alphaSigma_maxBound;
							isBadSigmaBound[k]= true;
							isSigmaBoundConstraintViolated= true;
							if(j==l) cout<<"INFO: Sigma max constraint violated for component no. "<<k+1<<":  Sigma(k)="<<Sigma<<"  Sigma_max(k)="<<Sigma_max<<"  Sigma_safe(k)="<<Sigma_safe<<"  alpha="<<alphaSigma_maxBound<<endl;
							else cout<<"INFO: Covariance max constraint violated for component no. "<<k+1<<":  Sigma(k)="<<Sigma<<"  Sigma_max(k)="<<Sigma_max<<"  Sigma_safe(k)="<<Sigma_safe<<"  alpha="<<alphaSigma_maxBound<<endl;

							if(alphaSigma_maxBound<0 || alphaSigma_maxBound>1){
								cerr<<"MSNMixtureFitter::EMConstrain(): ERROR: Invalid alpha value for covariance max bound constraint for component "<<k+1<<"-"<<k+2<<" (alpha="<<alphaSigma_maxBound<<") ...exit!"<<endl;
								exit(1);
							}
						}

					}//end loop dim
				}//end loop dim

				
				alphaSigmaBoundList.push_back(alphaSigmaBound);
			}//end loop components

		}//close if fUseCovarianceBoundConstraint


		//## COVARIANCE GROUP CONSTRAINT
		
		std::vector<double> alphaSigmaList;
		alphaSigmaList.clear();
		alphaSigmaList.resize(0);

		bool isBadSigma= false;
		bool isSigmaGroupConstraintViolated= false;

		if(fUseCovarianceConstraint){
			
			for(int j=0;j<fNDim;j++){
				double constraintSign= (*fSigmaConstraintSign)(j,0);
			
				for(int k=0;k<fNComponents-1;k++){
					
					double Sigma_k= fSigma[k](j,j);
					double Sigma_k_1= fSigma[k+1](j,j);
					double safeSigma_k= fSigma_safe[k](j,j);
					double safeSigma_k_1= fSigma_safe[k+1](j,j);
					double denomSigma= 1.-(Sigma_k_1-Sigma_k)/(safeSigma_k_1-safeSigma_k);
					double alphaSigma= 1./denomSigma;
		
					

					if(constraintSign*Sigma_k>constraintSign*Sigma_k_1){
						cout<<"MGMixtureFitter::EMConstrain(): INFO: Covariance constraint violated for component no. "<<k+1<<":  Sigma(k)="<<Sigma_k<<"  Sigma(k+1)="<<Sigma_k_1<<"  safeSigma(k)="<<safeSigma_k<<"  safeSigma(k+1)="<<safeSigma_k_1<<"  alpha="<<alphaSigma<<endl;

						if(alphaSigma<0 || alphaSigma>1){
							cerr<<"MGMixtureFitter::EMConstrain(): ERROR: Invalid alpha value for covariance group constraint for component "<<k+1<<"-"<<k+2<<" (alpha="<<alphaSigma<<") ...exit!"<<endl;
							exit(1);
						}

						isBadSigma= true;
						isSigmaGroupConstraintViolated= true;
						alphaSigmaList.push_back(alphaSigma);	
					}//close if Sigma constrain violated

				}//end loop components	
			}//end loop dim
		}//close if fUseCovarianceConstraint

	
		//## Computing min alpha
		cout<<"Computing min alpha "<<endl;
		
		double minAlpha_groupSigmaConstraint= 1;
		if(isSigmaGroupConstraintViolated){
			cout<<"MGMixtureFitter::EMConstrain(): INFO: alphaSigmaList(";
			for(unsigned int i=0;i<alphaSigmaList.size();i++){
				cout<<alphaSigmaList[i]<<",";
				if(alphaSigmaList[i]<minAlpha_groupSigmaConstraint) minAlpha_groupSigmaConstraint= alphaSigmaList[i];
			}	
			cout<<")"<<endl;
		}

		double minAlpha_boundSigmaConstraint[fNComponents];
		for(int k=0;k<fNComponents;k++){
			minAlpha_boundSigmaConstraint[k]= 1;
			if(isBadSigmaBound[k]){
				minAlpha_boundSigmaConstraint[k]= alphaSigmaBoundList[k].Min();	
			}
		}//end loop components

		

		if( isSigmaGroupConstraintViolated || isSigmaBoundConstraintViolated ){
			cout<<"MSNMixtureFitter::EMConstrain(): INFO: Updating covariance matrix..."<<endl;
			
			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateSigma= false;

			for(int k=0;k<fNComponents;k++){
				if(!isSigmaGroupConstraintViolated && !isBadSigmaBound[k]) continue;

				double alphaSigmaMin= min(minAlpha_groupSigmaConstraint,minAlpha_boundSigmaConstraint[k]);
				double alphaSigmaOpt= alphaSigmaMin/fConstraintAlphaScale;

				cout<<"MGMixtureFitter::EMConstrain(): INFO: Component no. "<<k+1<<": minAlpha_groupConstraint="<<minAlpha_groupSigmaConstraint<<"  minAlpha_boundConstraint="<<minAlpha_boundSigmaConstraint[k]<<"  alphaSigmaMin="<<alphaSigmaMin<<"  alphaSigmaOpt="<<alphaSigmaOpt<<endl;

				if(alphaSigmaOpt<fConstraintAlphaTolerance) {
					cerr<<"MGMixtureFitter::EMConstrain(): WARNING: Covariance for component "<<k+1<<" is stuck in constraint ("<<alphaSigmaOpt<<"<"<<fConstraintAlphaTolerance<<")... ...regenerate!"<<endl;
					hasToRegenerateSigma= true;
					hasToBeGenerated[k]= true;	
				}
			
					
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double Sigma= fSigma[k](j,l);
						double safeSigma= fSigma_safe[k](j,l);
						fSigma[k](j,l)= (1.-alphaSigmaOpt)*safeSigma + alphaSigmaOpt*Sigma;
					}//end loop dim
				}//end loop dim

			}//end loop components
		
			if(hasToRegenerateSigma && fUseRandomRegenerationAfterStuck){
				cout<<"MGMixtureFitter::EMConstrain(): INFO: Regenerating covariance pars..."<<endl;
				RandomizeSigmaParFromModel(hasToBeGenerated);	
			}
		
			//## Re-Calculate sigma starting from the updated covariance
			for(int k=0;k<fNComponents;k++){
				fSigmaEigen[k]= MathUtilities::GetEigenDecomposition(fSigma[k]);	
			}//end loop components

		}//close if (fUseCovarianceConstraint || fUseCovarianceBoundConstraint)
	

	}//close if fFixCovariancePar




	//## Set the current constrained set as "safe"
	for(int k=0;k<fNComponents;k++){
		fMu_safe[k]= fMu[k];
		fSigma_safe[k]= fSigma[k];
		fSigmaEigen_safe[k]= fSigmaEigen[k];	
		fP_safe[k]= fP[k];
		
		cout<<"MGMixtureFitter::EMConstrain(): Component "<<k+1<<" =="<<endl;
		cout<<"p= "<<fP[k]<<endl;
		cout<<"Mu= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(j,0)<<",";
		cout<<(fMu[k])(fNDim-1,0)<<")"<<endl; 
			
		cout<<"Sigma= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
	}//end loop components
	*/

	return 0;

}//close RunEM_ConstrainStep()



int MNMixtureClustering::InitParsToUser()
{
	cout<<"INFO: Initializing mixture parameters with user defaults..."<<endl;

	//Check fractions weights and set to starting values
	int nFractionPars= static_cast<int>(fOptions.P_start.size());
	if(nFractionPars!=fOptions.nComponents){
		cerr<<"ERROR: Number of components weights given as arg ("<<nFractionPars<<") is different from nComponents ("<<fOptions.nComponents<<")!"<<endl;
		return -1;
	}
	for(int k=0;k<fOptions.nComponents;k++){
		fP[k]= fOptions.P_start[k];
	}
	
	//Check component means
	cout<<"DEBUG: Print start means..."<<endl;
	for(size_t k=0;k<fOptions.Mu_start.size();k++){
		fOptions.Mu_start[k].Print();
	}

	int nMeanComponents= static_cast<int>(fOptions.Mu_start.size());
	if(nMeanComponents != fOptions.nComponents){
		cerr<<"ERROR: Number of user mean components ("<<nMeanComponents<<") is different from nComponents ("<<fOptions.nComponents<<")!"<<endl;
		return -1;
	}
	for(size_t k=0;k<(fOptions.Mu_start).size();k++){
		//int meanVectDim= (fOptions.Mu_start)[k].GetNrows();
		int meanVectDim= (fOptions.Mu_start)[k].GetNcols();
		if(meanVectDim!=fNDim){
			cerr<<"ERROR: User mean vector size for component no. "<<k+1<<" is not equal to nDim="<<fNDim<<"!"<<endl;
			return -1;
		}
		fMu[k]= (fOptions.Mu_start)[k];
	}//end loop components
		
	//Check sigmas
	int nSigmaComponents= static_cast<int>(fOptions.Sigma_start.size());
	if(nSigmaComponents != fOptions.nComponents){
		cerr<<"ERROR: Number of user sigma components ("<<nSigmaComponents<<") is different from nComponents ("<<fOptions.nComponents<<")!"<<endl;
		return -1;
	}
	for(size_t k=0;k<(fOptions.Sigma_start).size();k++){
		int nRows= (fOptions.Sigma_start)[k].GetNrows();
		int nCols= (fOptions.Sigma_start)[k].GetNcols();
		if(nRows!=fNDim || nCols!=fNDim){
			cerr<<"ERROR: User sigma matrix for component no. "<<k+1<<" has size different from nDim x nDim (nDim="<<fNDim<<")!"<<endl;
			return -1;
		}
		fSigma[k]= (fOptions.Sigma_start)[k];
	}//end loop components

	return 0;

}//close InitParsToUser()

int MNMixtureClustering::InitParsToKMeans()
{
	cout<<"INFO: Initializing mixture parameters with k-means..."<<endl;

	//## Perform Kmeans clustering on pre-completed data
	KMeansClustering kmeans;
	if(kmeans.RunKMedians(fRTableName_preImp,fOptions.nComponents)<0){
		cerr<<"ERROR: KMeans clustering run failed!"<<endl;
		return -1;
	}

	//## Retrieve data
	TMatrixD* clusterSizes= kmeans.GetClusterSizes();
	TMatrixD* clusterCenters= kmeans.GetClusterCenters();
	std::vector<TMatrixD*> clusterCovMatrixes= kmeans.GetClusterCovMatrixes();
	//cout<<"DEBUG: clusterSizes: "<<clusterSizes->GetNrows()<<" x "<<clusterSizes->GetNcols()<<endl;
	//cout<<"DEBUG: clusterCenters: "<<clusterCenters->GetNrows()<<" x "<<clusterCenters->GetNcols()<<endl;
	

	//## Assign cluster weights as fraction start
	cout<<"INFO: Assign the start fraction parameters according to kmeans cluster weights..."<<endl;
	for(int k=0;k<fOptions.nComponents;k++){
		double clusterSize= (*clusterSizes)(k,0);
		fP[k]= clusterSize/(double)(fN);
	}//end loop components

	//## Assign cluster covariance as covariance par start
	cout<<"INFO: Assign the start covariance parameters according to the kmeans cluster variances..."<<endl;
	for(int k=0;k<fOptions.nComponents;k++){		
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) {
				fSigma[k](j,l)= (*clusterCovMatrixes[k])(j,l);	
			}//close loop dim
		}//close loop dim	
	}//end loop components
	
	
	//## Assign cluster centers as mean start
	cout<<"INFO: Assign the start mean parameters according to the kmeans cluster centers ..."<<endl;
	for(int k=0;k<fOptions.nComponents;k++){
		for(int j=0;j<fNDim;j++){
			//fMu[k](j,0)= (*clusterCenters)(k,j);
			fMu[k](0,j)= (*clusterCenters)(k,j); 	
		}//end loop dim
	}//end loop components

	
	//## Print pars
	//cout<<"DEBUG: Printing par start values after kmeans..."<<endl;
	//PrintPars();

	return 0;

}//close InitParsToKMeans()


void MNMixtureClustering::RandomizePars()
{
	cout<<"Randomize initial fit parameters..."<<endl;

	//## Generate random Sigma
	if(fOptions.randomizeStartCovariancePars) {
		cout<<"INFO: Randomizing mixture covariance pars..."<<endl;
		RandomizeSigmaPars(); 
	}

	//## Generate random means
	if(fOptions.randomizeStartMeanPars) {
		cout<<"INFO: Randomizing mixture mean pars..."<<endl;
		RandomizeMeanPars();
	}

	//## Print init parameters
	PrintPars();

}//close RandomizePars()

void MNMixtureClustering::RandomizeMeanPars()
{	
	//## Generate mu	
	for(int j=0;j<fNDim;j++){
		double parValue[fOptions.nComponents];
		double meanSigmaGeneration= 3.*sqrt((*fCovarianceData)(j,j));	
		//double meanSigmaGeneration= 3.*sqrt((*fSigma[k])(j,j));	
		//double minBoundary= (*fMeanData)(j,0)-meanSigmaGeneration;
		//double maxBoundary= (*fMeanData)(j,0)+meanSigmaGeneration;
		double minBoundary= (*fMeanData)(0,j)-meanSigmaGeneration;
		double maxBoundary= (*fMeanData)(0,j)+meanSigmaGeneration;

		for(int k=0;k<fOptions.nComponents;k++){
			parValue[k]= gRandom->Uniform(minBoundary,maxBoundary);
			//fMu[k](j,0)= parValue[k];
			fMu[k](0,j)= parValue[k];
		}//end loop components
	}//end loop dim
	
}//close RandomizeMeanPars()


void MNMixtureClustering::RandomizeSigmaPars()
{
	//## Generate Sigma
	TMatrixD randSigma(fNDim,fNDim);
	randSigma.Zero();
	TMatrixD S(fNDim,fNDim);
	S= *fCovarianceData;
	TMatrixD DiagS(fNDim,fNDim);
	DiagS= *fCovarianceDiagData;

	for(int k=0;k<fOptions.nComponents;k++){
		
		bool isSigmaPosDef= false;	
		while(!isSigmaPosDef){
			double rand= gRandom->Uniform(0,1);
	
			randSigma.Zero();
			randSigma= S + (rand-1)*DiagS;

			fSigma[k]= randSigma;

			double Det= fSigma[k].Determinant();
			if(Det>0) {	
				isSigmaPosDef= true;	
				//fInitRandNumbers[k]= rand;
			}
		}//end loop while
	}//end loop components
	
}//close RandomizeSigmaPars()


int MNMixtureClustering::CheckCovariance()
{
	//## Loop over components and check covariance matrix:
	//##   - must be symmetric
	//##   - must be positive def	
	//##   - set diagonal (if enabled)

	for(int k=0;k<fOptions.nComponents;k++)
	{
		//## Force diagonal covariance?
		if(fOptions.forceDiagonalCovariance){
			Util::MakeDiagonalMatrix(fSigma[k]);
		}

		//## Force covariance to be symmetric and pos def
		//## Approximate to the nearest covariance matrix
		if(Util::MakeSymmPosDefCovarianceMatrix(fSigma[k])<0){
			cerr<<"ERROR: Failed to correct covariance matrix for component no. "<<k+1<<"!"<<endl;
			return -1;	
		}

		//## Check if determinant is ok
		if(fSigmaDet[k]<=1.e-16){
			cout<<"WARN: Covariance matrix for component "<<k+1<<" is not positive def (det="<<fSigmaDet[k]<<") ... setting covariance to safe matrix!"<<endl;
			fSigma[k]= fSigma_safe[k];
		}


		//## Update inverse & determinant	
		fSigmaInv[k]= TMatrixD(TMatrixD::kInverted,fSigma[k]);
		fSigmaDet[k]= fSigma[k].Determinant();	
		
		//## Update eigenvalues & eigenvectors	
		if(Util::ComputeSymMatrixEigenvalues(fSigmaEigen[k],fSigmaEigenvect[k],fSigma[k])<0){
			cerr<<"ERROR: Failed to compute matrix eigenvalues/eigenvectors!"<<endl;
			return -1;
		}
		TMatrixD sigmaEigenvectInv= TMatrixD(TMatrixD::kInverted,fSigmaEigenvect[k]);
		double sigmaEigenvectDet= sigmaEigenvectInv.Determinant();
		if(sigmaEigenvectDet<=0) {
			cerr<<"WARN: Covariance matrix eigenvectors inversion failed for component "<<k+1<<" (det="<<sigmaEigenvectDet<<")!"<<endl;
		}

		//## Set the current covariance as "safe"
		fSigma_safe[k]= fSigma[k];
		fSigmaEigen_safe[k]= fSigmaEigen[k];
	
	}//end loop components

	return 0;

}//close CheckCovariance()



void MNMixtureClustering::PrintPars()
{
	for(int k=0;k<fOptions.nComponents;k++){
		cout<<"== Component "<<k+1<<" =="<<endl;
		cout<<"p= "<<fP[k]<<endl;
		cout<<"Mu= (";
		//for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(j,0)<<",";
		//cout<<(fMu[k])(fNDim-1,0)<<")"<<endl; 
		for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(0,j)<<",";
		cout<<(fMu[k])(0,fNDim-1)<<")"<<endl; 
				
		cout<<"Sigma= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
	}//end loop components

}//close PrintPars()


void MNMixtureClustering::PrintPars2()
{
	TString meanPrintMsg= "";
	TString sigmaPrintMsg= "";
	TString fractionPrintMsg= "p=(";

	for(int k=0;k<fOptions.nComponents;k++)
	{		
		fractionPrintMsg+= Form("%1.2f",fP[k]);
		if(k==fOptions.nComponents-1) fractionPrintMsg+= TString(")");
		else fractionPrintMsg+= TString(",");

		meanPrintMsg+= Form("Mu%d=(",k+1);
		//for(int j=0;j<fNDim-1;j++) meanPrintMsg+= Form("%1.2f,",(fMu[k])(j,0));
		//meanPrintMsg+= Form("%1.2f) ",(fMu[k])(fNDim-1,0));
		for(int j=0;j<fNDim-1;j++) meanPrintMsg+= Form("%1.2f,",(fMu[k])(0,j));
		meanPrintMsg+= Form("%1.2f) ",(fMu[k])(0,fNDim-1));

		sigmaPrintMsg+= Form("Sigma%d=(",k+1);
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) {
				if(j==fNDim-1 && l==fNDim-1) sigmaPrintMsg+= Form("%1.2f) ",(fSigma[k])(j,l));
				else sigmaPrintMsg+= Form("%1.2f,",(fSigma[k])(j,l));
			}
		}
	}//end loop components
	
	cout<<fractionPrintMsg<<endl;	
	cout<<meanPrintMsg<<endl;
	cout<<sigmaPrintMsg<<endl;
	
}//close PrintPars2()


}//close namespace

