/**
* @file Util.h
* @class Util
* @brief Util functions for 
*
* Useful math functions
* @author S. Riggi
* @date 17/09/2012
*/



#ifndef _UTIL_h
#define _UTIL_h 1


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

class Util : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Util();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~Util();

		static RInside fR;
		
	
	public:

		/**
		* \brief Create artifically random missing data in matrix and return data with missing data
		*/
		static TMatrixD* MakeRandomMissingData(TMatrixD* dataMatrix,double missingDataFraction);
		

		static int DumpMatrixToAsciiFile(TMatrixD* dataMatrix,std::string filename);

		/**
		* \brief Order vectors and get ordering index
		*/
		template<class T> struct index_cmp{

  		index_cmp(const T arr) : arr(arr) {}
  		bool operator()(const size_t a, const size_t b) const
 			{
    		return arr[a] > arr[b];
  		}
  		const T arr;
		};

		template< class T >
			static void reorder(std::vector<T> & unordered,std::vector<size_t> const & index_map,std::vector<T> & ordered){
  			// copy for the reorder according to index_map, because unsorted may also be
  			// sorted
  			std::vector<T> copy = unordered;
  			ordered.resize(index_map.size());
  			for(int i = 0; i<index_map.size();i++)
					ordered[i] = copy[index_map[i]];
			}

		template <class T>
			static void sort(std::vector<T> & unsorted,std::vector<T> & sorted,std::vector<size_t> & index_map){
  			// Original unsorted index map
  			index_map.resize(unsorted.size());
 				for(size_t i=0;i<unsorted.size();i++)
					index_map[i] = i;
  
  			// Sort the index map, using unsorted for comparison
  			std::sort(index_map.begin(),index_map.end(),index_cmp<std::vector<T>& >(unsorted));
  			sorted.resize(unsorted.size());
  			reorder(unsorted,index_map,sorted);
			}

	private:
	
	

		friend class MeanImputation;
		friend class MultipleImputation;
		friend class DataReader;

	ClassDef(Util,1)

};//close class

}//close namespace

#endif




