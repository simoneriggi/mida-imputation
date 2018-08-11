/**
* @file MeanImputation.h
* @class MeanImputation
* @brief MeanImputation
*
* Perform mean imputation on a data sample
* @author S. Riggi
* @date 30/07/2013
*/



#ifndef _MEAN_IMPUTATION_h
#define _MEAN_IMPUTATION_h 1


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

class MeanImputation{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MeanImputation();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MeanImputation();

	
	public:
		/**
		* \brief Run imputation
		*/
		static TMatrixD* RunImputation(std::string filename,std::string delimiter="");
		/**
		* \brief Run imputation
		*/
		static TMatrixD* RunImputation(TMatrixD* dataMatrix);

		
	private:
	
	private:
	

	ClassDef(MeanImputation,1)


};//close class

}//close namespace 

#endif
