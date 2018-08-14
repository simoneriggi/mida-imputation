/**
* @file ListwiseDeletion.h
* @class ListwiseDeletion
* @brief ListwiseDeletion
*
* Perform list-wise deletion on a data sample
* @author S. Riggi
* @date 30/07/2013
*/



#ifndef _LISTWISE_DELETION_h
#define _LISTWISE_DELETION_h 1


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


class ListwiseDeletion : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    ListwiseDeletion();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~ListwiseDeletion();

	
	public:
	
		/**
		* \brief Run imputation
		*/
		static TMatrixD* RunImputation(std::string filename,std::string delimiter="");
		/**
		* \brief Run imputation
		*/
		static TMatrixD* RunImputation(TMatrixD* dataMatrix);
		/**
		* \brief Run imputation from loaded R table 
		*/
		static TMatrixD* RunImputationFromRTable(std::string RTableName,std::string RTableName_withoutMiss= "dataMatrix_imputed");
		/**
		* \brief Run imputation from loaded R table and return results in a R table
		*/
		static int RunImputationInR(std::string RTableName,std::string RTableName_withoutMiss);

	private:


		ClassDef(ListwiseDeletion,1)

};//close class

}//close namespace

#endif

