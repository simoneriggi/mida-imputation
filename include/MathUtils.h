/**
* @file MathUtils.h
* @class MathUtils
* @brief Mathematical utility functions
*
* Useful math functions
* @author S. Riggi
* @date 17/09/2012
*/

#ifndef _MATH_UTILS_h
#define _MATH_UTILS_h 1


#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TF1.h>
#include <TGraph.h>
#include <TVector3.h>

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


class MathUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MathUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MathUtils();

		

	public:
		/**
		* \brief Multivariate gaussian function
		*/
		static double MGDensityFcn(TMatrixD& y,TMatrixD& mu,TMatrixD& SigmaInv,double SigmaDet);

		/**
		* \brief Compute posterion probability 
		*/
		static double tauGaus(TMatrixD& y,TMatrixD& mu,TMatrixD& SigmaInv,double SigmaDet,double pi);
		
	private:
	
		friend class MNMixtureClustering;

	ClassDef(MathUtils,1)

};//close class

}//close namespace 

#endif

