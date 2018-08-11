/**
* @file DataReader.h
* @class DataReader
* @brief DataReader
*
* Read ascii data with missing values and export to matrix
* @author S. Riggi
* @date 30/07/2013
*/

#ifndef _DATA_READER_h
#define _DATA_READER_h 1

#include <TMatrixD.h>

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

class DataReader : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    DataReader();
		/**
		* \brief Class destructor: free allocated memory
		*/
   virtual ~DataReader();

	
	public:
		/**
		* \brief Read ascii data and export to TMatrixD
		*/
		static TMatrixD* ReadAscii(std::string filename,std::string delimiter="");

		/**
		* \brief Import data matrix in a R table
		*/
		static int ReadAsciiInR(std::string filename,std::string delimiter="",std::string importedDataName="data");


		/**
		* \brief Read ascii data 
		*/
		static int ReadAscii(std::vector< std::vector<double> >& data,long int& ndim, std::string filename,std::string delimiter="");

	ClassDef(DataReader,1)

};//close class DataReader

}//close namespace

#endif

