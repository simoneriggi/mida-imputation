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

#include <Logger.h>
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
		static double MNDensityFcn(TMatrixD& y,TMatrixD& mu,TMatrixD& SigmaInv,double SigmaDet);

		/**
		* \brief Compute posterion probability 
		*/
		static double tauGaus(TMatrixD& y,TMatrixD& mu,TMatrixD& SigmaInv,double SigmaDet,double pi);
		
	private:
	
		friend class MNMixtureClustering;

	ClassDef(MathUtils,1)

};//close class

//===============================
//==  MULTIVARIATE GAUS FCN
//===============================
class MNPDF : public TObject {
 
	public:
		TMatrixD X;
		TMatrixD Mu;
		TMatrixD CovMat; 
   
	public:
   	MNPDF(int dim) : 
      X(TMatrixD(1,dim)), 
      Mu(TMatrixD(1,dim)), 
      CovMat(TMatrixD(dim,dim))
   	{}

		virtual ~MNPDF(){}
   
	public:

		double operator() (double *x, double *p) 
		{ 
			int dim= X.GetNcols();
			int parCounter= 0;

			//Amplitude
			double A= p[parCounter];
			parCounter++;			

			//Fill data & mean vectors
			for (int i=0;i<dim;++i) { 
				X(0,i)= x[i];
				Mu(0,i)= p[parCounter]; 
				parCounter++; 
			} 

			//Fill covariance matrix
			for (int i=0; i<dim; ++i) { 
      	for (int j=i; j<dim; ++j) { 
					CovMat(i,j)= p[parCounter];
					CovMat(j,i)= p[parCounter];
					parCounter++;
				}
			}

			//Print 
			//X.Print();
			//Mu.Print();
      //CovMat.Print();

			//Compute inverse
			TMatrixD CovMatInv= TMatrixD(TMatrixD::kInverted,CovMat);

			double det = CovMat.Determinant(); 
      if (det <= 0) {
      	ERROR_LOG("Determinant is <=0 (det="<<det<<"), returning 0!");
      	return 0;  
      }

			//Compute function value
			double fval= A*MathUtils::MNDensityFcn(X,Mu,CovMatInv,det);

      return fval;

	}//close function

	ClassDef(MNPDF,1)

};//close class

//===================================
//==  MULTIVARIATE GAUS MIXTURE FCN
//===================================
class MNMixturePDF : public TObject {
 
	public:
		int nDim;
		int nComponents;
		int nComponentPars;
   
	public:
   	MNMixturePDF(int dim,int k) : 
      nDim(dim), nComponents(k)
   	{
			//Compute number of component pars
			//amplitude + nDim means + nDim*(nDim+1)/2 cov elements
			nComponentPars= 1 + nDim + nDim*(nDim+1)/2;
		}//close constructor
   
		virtual ~MNMixturePDF(){}

	public:

		double operator() (double *x, double *p) 
		{ 
			double fval= 0.;
			double componentPars[nComponentPars];
			int parCounter= 0;
			for(int k=0;k<nComponents;k++){
				for(int j=0;j<nComponentPars;j++){
					componentPars[j]= p[parCounter];
					parCounter++;	
				}//end loop component pars

				//cout<<"componentPars {";
				//for(int j=0;j<nComponentPars;j++) cout<<componentPars[j]<<",";
				//cout<<"}"<<endl;
				
				MNPDF fcn(nDim);
				double fval_component= fcn(x,componentPars); 
				fval+= fval_component;
			}//end loop nComponents

      return fval;

	}//close function

	ClassDef(MNMixturePDF,1)

};//close class


}//close namespace 

#endif

