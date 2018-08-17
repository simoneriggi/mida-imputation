/**
* @file DataGenerator.h
* @class DataGenerator
* @brief Data generator for multivariate normal/skew-normal/skew-t data
*
* Random generator
* @author S. Riggi
* @date 17/09/2012
*/

#ifndef _DATA_GENERATOR_h
#define _DATA_GENERATOR_h 1


#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>


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

class DataGenerator : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    DataGenerator();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~DataGenerator();

		
	public:

		/**
		* \brief Generate data sample from multivariate normal 
		*/
		static TMatrixD* GenMNSample (
			int N,
			std::vector<double>& minRange,
			std::vector<double>& maxRange,
			TMatrixD& mean,
			TMatrixD& covMatrix,
			bool addMissingData=false,
			double missDataFraction=0
		);

		/**
		* \brief Generate data sample from multivariate normal 
		*/
		static int GenMNSample (
			TMatrixD* dataMatrix,
			int N,
			std::vector<double>& minRange,
			std::vector<double>& maxRange,
			TMatrixD& mean,
			TMatrixD& covMatrix,
			int dataFillOffset=0,
			bool addMissingData=false,
			double missDataFraction=0
		);

		/**
		* \brief Generate data sample from multivariate normal mixture
		*/
		static TMatrixD* GenMNMixtureSample(	
			int N,
			std::vector<double>& minRange,
			std::vector<double>& maxRange,
			std::vector<double>& weights,
			std::vector<TMatrixD>& means,
			std::vector<TMatrixD>& covMatrixes,
			bool addMissingData=false,
			double missDataFraction=0
		);

		

	private:
		
		
	ClassDef(DataGenerator,1)


};//close class

}//close namespace

#endif


