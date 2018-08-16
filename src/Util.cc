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

int Util::LoadRLibraries(std::vector<std::string> libraryNames)
{
	//Check libraries
	if(libraryNames.empty()){
		cerr<<"WARN: Empty library names!"<<endl;
		return -1;
	}

	//Load libraries
	for(size_t i=0;i<libraryNames.size();i++){
		std::string RCmd= Form("library('%s');",libraryNames[i].c_str());
		try{
			Util::fR.parseEval(RCmd);	
		}
		catch(...){
			cerr<<"ERROR: Failed to load library "<<libraryNames[i]<<"!"<<endl;
			return -1;
		}		
	}//end loop libraries

	return 0;

}//close LoadRLibraries()

int Util::ClearRData()
{
	//## Clear R environment
	cout<<"INFO: Clearing R environment..."<<endl;
	try{
		fR.parseEvalQ("rm(list = ls(all = TRUE));");
	}
	catch(...){
		cerr<<"ERROR: Failed to clear R data!"<<endl;
		return -1;
	}

	return 0;

}//close ClearRData()


int Util::MakeSymmetricMatrix(TMatrixD& C)
{
	//Check if already symmetric
	if(C.IsSymmetric()){
		cout<<"INFO: Matrix is already symmetric, nothing to be done..."<<endl; 
		return 0;
	}

	//Force matrix symmetric
	TMatrixD C_t(TMatrixD::kTransposed,C);
	TMatrixD C_sym = 0.5*(C + C_t);
	C= C_sym;

	return 0;

}//close MakeSymmetricMatrix()


int Util::ComputeSymMatrixEigenvalues(TMatrixD& eigenVals,TMatrixD& eigenVects,const TMatrixD& C)
{
	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();

	//Resize eigen matrix
	eigenVals.ResizeTo(1,nCols);
	eigenVects.ResizeTo(nRows,nCols);

	//Check for square matrix
	if(nCols!=nRows){
		cerr<<"ERROR: Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)"<<endl;
		return -1;
	}
	
	//Check symmetric
	if(!C.IsSymmetric()){
		cerr<<"ERROR: Input matrix is not symmetric!"<<endl;
		return -1;
	}

	//Compute eigenvalues & eigenvect
	TMatrixDEigen matrixDecomposition(C);
	TVectorD eigenVals_re= matrixDecomposition.GetEigenValuesRe();
	eigenVects= matrixDecomposition.GetEigenVectors();
	for(int i=0;i<eigenVals_re.GetNrows();i++) eigenVals(0,i)= eigenVals_re(i);

	return 0;

}//close ComputeMatrixEigenvalues()

int Util::ComputeMatrixEigenvalues(TMatrixD& eigenVals,TMatrixD& eigenVects,TMatrixD& C,bool forceSymmetric)
{
	//## NB: If matrix is not simmetric then the eigenvalue matrix D is block
	//       diagonal with the real eigenvalues in 1-by-1 blocks and any complex
	//       eigenvalues, a + i*b, in 2-by-2 blocks, [a, b; -b, a]
	//       If matrix is simmetric the eigenvalue matrix is diagonal with real eigenvalue on the diagonal

	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();

	//Resize eigen matrix
	eigenVals.ResizeTo(nRows,nCols);
	eigenVects.ResizeTo(nRows,nCols);

	//Check for square matrix
	if(nCols!=nRows){
		cerr<<"ERROR: Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)"<<endl;
		return -1;
	}
	
	//Make symmetric
	if(forceSymmetric && MakeSymmetricMatrix(C)){
		cerr<<"ERROR: Failed to force simmetric matrix!"<<endl;
		return -1;
	}

	//Compute eigenvalues & eigenvect
	TMatrixDEigen matrixDecomposition(C);
	eigenVals= matrixDecomposition.GetEigenValues();
	eigenVects= matrixDecomposition.GetEigenVectors();

	return 0;

}//close ComputeMatrixEigenvalues()

int Util::ComputeMatrixEigenvaluesInR(TMatrixD& eigenVals,TMatrixD& eigenVects,TMatrixD& C,bool forceSymmetric)
{
	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();
	
	//Resize eigen matrix
	eigenVals.ResizeTo(nCols,1);
	eigenVects.ResizeTo(nRows,nCols);

	//Check for square matrix
	if(nCols!=nRows){
		cerr<<"ERROR: Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)"<<endl;
		return -1;
	}

	//Import matrix in R
	std::string matrixRName= "sigma";
	if(ImportMatrixInR(&C,matrixRName)<0){
		cerr<<"ERROR: Failed to import matrix in R!"<<endl;
		return -1;
	}

	//Force symmetric?
	if(forceSymmetric){
		try{
			fR.parseEvalQ(Form("forceSymmetric(%s);",matrixRName.c_str()));
		}
		catch(...){
			cerr<<"ERROR; Failed to force symmetric matrix in R!"<<endl;
			return -1;
		}
	}//close if

	//Compute eigenvalues & vectors
	std::stringstream ss;
	ss<<"eigenRes <- eigen("<<matrixRName<<");";
	std::string RCmd= ss.str();
	try{
		fR.parseEvalQ(RCmd);
		fR.parseEvalQ(Form("rm(%s);",matrixRName.c_str()));//remove tmp sigma
	}
	catch(...){
		cerr<<"ERROR: Failed to compute eigenvalues & eigenvectors of given matrix!"<<endl;
		return -1;
	}

	//Retrieve results
	try{
		Rcpp::NumericVector sigmaEigenValues= fR.parseEval("eigenRes$values;");
		Rcpp::NumericMatrix sigmaEigenVectors= fR.parseEval("eigenRes$vectors;");

		for (int l=0; l<nCols; l++) {
			eigenVals(0,l)= sigmaEigenValues(l);
			for (int j=0; j<nCols; j++) {
				double w= sigmaEigenVectors(l,j);
				eigenVects(l,j)= w;
			}
		}
		
		fR.parseEvalQ("rm(eigenRes);");

	}//close try
	catch(...){
		cerr<<"ERROR: Failed to retrieve eigenvalues/eigenvector matrix!"<<endl;
		return -1;
	}

	return 0;
	
}//close ComputeMatrixEigenvalues()


void Util::MakeDiagonalMatrix(TMatrixD& C)
{
	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();
	
	for(int j=0;j<nRows;j++){
		for(int l=0;l<nCols;l++){
			if(j==l) continue;
			C(j,l)= 0.;
		}//end loop dim
	}//end loop dim

}//close MakeDiagonalMatrix()


int Util::MakeSymmPosDefCovarianceMatrix(TMatrixD& covMatrix)
{
	long int nCols= covMatrix.GetNcols();
	long int nRows= covMatrix.GetNrows();
	
	//Check for square matrix
	if(nCols!=nRows){
		cerr<<"ERROR: Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)"<<endl;
		return -1;
	}

	//Import matrix in R
	std::string covMatrixRName= "sigma";
	if(ImportMatrixInR(&covMatrix,covMatrixRName)<0){
		cerr<<"ERROR: Failed to import covariance matrix in R!"<<endl;
		return -1;
	}

	//## Force covariance to be symmetric and pos def
	//## Approximate to the nearest covariance matrix
	std::stringstream ss;
	ss<<"res <- nearPD("<<covMatrixRName<<", corr=FALSE, do2eigen=FALSE, ensureSymmetry= TRUE);";
	std::string RCmd= ss.str();

	try{
		fR.parseEvalQ(RCmd);

		//Remove tmp sigma
		fR.parseEvalQ(Form("rm(%s);",covMatrixRName.c_str()));
	}
	catch(...){
		cerr<<"ERROR: Failed to approximate covariance to nearest symm & pos-def matrix!"<<endl;	
		return -1;
	}
	
	//## Get corrected matrix and re-assign to Sigma
	try{
		Rcpp::NumericMatrix covMatrix_corr= fR.parseEval("as.matrix(res$mat);");
		for (int l=0; l<nCols; l++) {
			for (int j=0; j<nCols; j++) {
				double w= covMatrix_corr(l,j);
				covMatrix(l,j)= w;
			}
		}
	}//close try
	catch(...){
		cerr<<"ERROR: Failed to retrieve corrected cov matrix and update given matrix!"<<endl;
		return -1;
	}

	return 0;

}//close MakeSymmPosDefCovarianceMatrix()

TMatrixD* Util::GetDiagonalMatrix(TMatrixD* dataMatrix)
{
	//Check data matrix
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to matrix given!"<<endl;	
		return nullptr;
	}
	long int nCols= dataMatrix->GetNcols();
	long int nRows= dataMatrix->GetNrows();
	
	//Check for square matrix
	if(nCols!=nRows){
		cerr<<"ERROR: Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a simmetric matrix!)"<<endl;
		return nullptr;
	}

	//Create matrix
	TMatrixD* diagMatrix= new TMatrixD(nRows,nRows);
	diagMatrix->Zero();

	//Fill diagonal values
	for(long int i=0;i<nRows;i++){
		(*diagMatrix)(i,i)= (*dataMatrix)(i,i);
	}

	return diagMatrix;

}//close GetDiagonalMatrix()


TMatrixD* Util::ComputeRTableColMeans(std::string RTable,std::string colMeansRName)
{
	//Check table name
	if(RTable==""){
		cerr<<"ERROR: Empty R table name!"<<endl;	
		return nullptr;
	}

	//Compute col means matrix	
	std::stringstream ss;
	ss<<colMeansRName<<" <- colMeans("<<RTable<<")";
	std::string RCmd= ss.str(); 
	try{
		Util::fR.parseEval(RCmd.c_str());
	}
	catch(...){
		cerr<<"ERROR: Failed to compute data column means in R!"<<endl;
		return nullptr;
	}

	//Convert data to TMatrixD
	TMatrixD* colMeans= Util::ConvertRVectToROOTMatrix(colMeansRName);
	if(!colMeans){
		cerr<<"ERROR: Failed to convert R vector to ROOT!"<<endl;
		return nullptr;
	}

	return colMeans;

}//close ComputeRTableColMeans()


TMatrixD* Util::ComputeCovarianceMatrixFromRTable(std::string RTable,std::string covMatrixRName)
{
	//Check table name
	if(RTable==""){
		cerr<<"ERROR: Empty R table name!"<<endl;	
		return nullptr;
	}

	//Compute covariance matrix	
	std::stringstream ss;
	ss<<covMatrixRName<<" <- cov("<<RTable<<")";
	std::string RCmd= ss.str(); 

	try{
		Util::fR.parseEval(RCmd.c_str());
	}
	catch(...){
		cerr<<"ERROR: Failed to compute data covariance matrix in R!"<<endl;
		return nullptr;
	}

	//Convert data to TMatrixD
	TMatrixD* covMatrix= Util::ConvertRTableToROOTMatrix(covMatrixRName);
	if(!covMatrix){
		cerr<<"ERROR: Failed to convert cov matrix from R to ROOT!"<<endl;
		return nullptr;
	}

	return covMatrix;

}//close ComputeCovarianceMatrixFromRTable()


TMatrixD* Util::ConvertRVectToROOTMatrix(std::string RVect)
{
	//Check table name
	if(RVect==""){
		cerr<<"ERROR: Empty R vector name!"<<endl;	
		return nullptr;
	}

	//Store R table to NumericVector and then fill TMatrixD (double copy...not efficient!!!)
	TMatrixD* dataMatrix_ROOT= 0;
	try{
		Rcpp::NumericVector dataVect= Util::fR.parseEval(RVect.c_str());
		long int N= Util::fR.parseEval(Form("length(%s)",RVect.c_str()));
		//dataMatrix_ROOT= new TMatrixD(N,1);
		dataMatrix_ROOT= new TMatrixD(1,N);

		for(long int i=0;i<N;i++){
			//(*dataMatrix_ROOT)(i,0)= dataVect(i);
			(*dataMatrix_ROOT)(0,i)= dataVect(i);
		}
	}//close try block
	catch(...){
		cerr<<"ERROR: Failed to retrieve data table and relative size with imputed values in R!"<<endl;
		return nullptr;
	}

	return dataMatrix_ROOT;

}//close ConvertRVectToROOTMatrix()


TMatrixD* Util::ConvertRTableToROOTMatrix(std::string RTable)
{
	//Check table name
	if(RTable==""){
		cerr<<"ERROR: Empty R table name!"<<endl;	
		return nullptr;
	}

	//Store R table to Numeric matrix and then fill TMatrixD (double copy...not efficient!!!)
	TMatrixD* dataMatrix_ROOT= 0;
	try{
		Rcpp::NumericMatrix dataMatrix= Util::fR.parseEval(RTable.c_str());
		
		//Util::fR.parseEvalQ("N <- nrow(dataMatrix);");
		//Util::fR.parseEvalQ("NDim <- ncol(dataMatrix);");
		//long int N= Util::fR.parseEval("N");
		//long int NDim= Util::fR.parseEval("NDim");
		long int N= Util::fR.parseEval(Form("nrow(%s)",RTable.c_str()));
		long int NDim= Util::fR.parseEval(Form("ncol(%s)",RTable.c_str()));
		
		dataMatrix_ROOT= new TMatrixD(N,NDim);
		for(long int i=0;i<N;i++){
			for(long int j=0;j<NDim;j++){
				(*dataMatrix_ROOT)(i,j)= dataMatrix(i,j);
			}
		}
	}//close try block
	catch(...){
		cerr<<"ERROR: Failed to retrieve data table and relative size with imputed values in R!"<<endl;
		return nullptr;
	}

	return dataMatrix_ROOT;

}//close ConvertRTableToROOTMatrix()

int Util::ImportMatrixInR(TMatrixD* dataMatrix,std::string dataname)
{
	//## Comvert matrix to R
	Rcpp::NumericMatrix* matrix_r= ConvertROOTMatrixToRMatrix(dataMatrix);
	if(!matrix_r){
		cerr<<"ERROR: Failed to convert ROOT matrix to R!"<<endl;
		return -1;
	}

	//## Import in R prompt
	try{
		fR[dataname.c_str()]= *matrix_r;
	}
	catch(...){
		cerr<<"ERROR: Failed to import RNumeric matrix in R prompt!"<<endl;
		return -1;
	}

	//Test 
	//std::stringstream ss;
	//ss<<"print("<<dataname<<")";
	//std::string RCmd= ss.str();
	//Util::fR.parseEvalQ(RCmd.c_str());
	
	return 0;

}//close ImportMatrixInR()

Rcpp::NumericMatrix* Util::ConvertROOTMatrixToRMatrix(TMatrixD* dataMatrix)
{
	//## Check data
	if(!dataMatrix){
		cerr<<"ERROR: Null ptr to data matrix given!"<<endl;
		return nullptr;
	}

	//## Create NumericMatrix
	long int nDim= dataMatrix->GetNcols();
	long int N= dataMatrix->GetNrows();
	Rcpp::NumericMatrix* matrix_r= new Rcpp::NumericMatrix(N,nDim);
	for(long int i=0;i<N;i++){
		for(long int j=0;j<nDim;j++){
			(*matrix_r)(i,j)= (*dataMatrix)(i,j);
		}//end loop dim
	}//end loop events

	return matrix_r;

}//close ConvertROOTMatrixToRMatrix()


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
