/**
* @file MultipleImputation.h
* @class MultipleImputation
* @brief MultipleImputation
*
* Perform multiple imputation on a data sample
* @author S. Riggi
* @date 30/07/2013
*/



#ifndef _MULTIPLE_IMPUTATION_h
#define _MULTIPLE_IMPUTATION_h 1


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


class MultipleImputation : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MultipleImputation();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MultipleImputation();

	
	public:
		/**
		* \brief Run imputation
		*/
		static TMatrixD* RunImputation(std::string filename,int nRuns=5,int nRepeatedRuns=100,std::string delimiter="");
		/**
		* \brief Run imputation
		*/
		static TMatrixD* RunImputation(TMatrixD* dataMatrix,int nRuns=5,int nRepeatedRuns=100);


	private:

	
		ClassDef(MultipleImputation,1)

};//close class

}//close namespace 

#endif

