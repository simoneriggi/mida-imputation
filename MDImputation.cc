#include <MeanImputation.h>
#include <ListwiseDeletion.h>
#include <MultipleImputation.h>
#include <Util.h>

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
	eMN= 4,
	eMSN= 5
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
  {(char*)0, (int)0, (int*)0, (int)0}
};

//Options
int imputationMethod= 1;
std::string inputFileName= "";
std::string outputFileName= "Output.dat";
long int nRuns= 5;
long int nRepeatedRuns= 100;

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

	while((c = getopt_long(argc, argv, "hi:o::m::r::R::",options_tab, &option_index)) != -1) {
    
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
		inputedDataMatrix= MeanImputation::RunImputation(inputFileName);
	}
	else if(imputationMethod==eLD){
		//inputedDataMatrix= ListwiseDeletion::RunImputation(inputFileName,{'\n'});
		//inputedDataMatrix= ListwiseDeletion::RunImputation(inputFileName,{' '});
		inputedDataMatrix= ListwiseDeletion::RunImputation(inputFileName," ");
	}
	else if(imputationMethod==eMI){
		inputedDataMatrix= MultipleImputation::RunImputation(inputFileName,nRuns,nRepeatedRuns);
	}
	else{
		cerr<<"ERROR: Invalid imputation methodselected (method="<<imputationMethod<<") ...exit!"<<endl;
		exit(1);
	}

	//## Check output
	if(!inputedDataMatrix){
		cerr<<"ERROR: Failed to compute imputed data!"<<endl;
		exit(1);
	}

	//## Save data to file
	if(Util::DumpMatrixToAsciiFile(inputedDataMatrix,outputFileName)<0){
		cerr<<"ERROR: Failed to dump imputed data to file!"<<endl;
		exit(1);
	}

	return 0;
	
}//close macro


