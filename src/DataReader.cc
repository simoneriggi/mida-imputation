/**
* @file DataReader.cc
* @class DataReader
* @brief DataReader
*
* Read ascii data with missing values and export to matrix
* @author S. Riggi
* @date 30/07/2013
*/

#include <DataReader.h>
#include <Util.h>

#include <TMatrixD.h>
#include <TMath.h>

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
#include <climits>

using namespace std;

ClassImp(MDImputation_ns::DataReader)

namespace MDImputation_ns {

DataReader::DataReader()
{

}

DataReader::~DataReader()
{

}

int DataReader::ReadAsciiInR(std::string filename,std::string delimiter,std::string importedDataName)
{
	//## Create R command for importing data
	std::stringstream RCmdStr;
	RCmdStr<< importedDataName<<" <- as.matrix(read.table(";
	RCmdStr<< "\""<< filename <<"\"";
	RCmdStr<< ", sep=\""<<delimiter<<"\"";
	RCmdStr<< "))";

	std::string RCmd= RCmdStr.str();
	cout<<"INFO: Running command in R: "<<RCmd<<endl;

	//## Read data matrix
	try {
		Util::fR.parseEvalQ(RCmd);
	}
	catch(...){
		cerr<<"ERROR: Failed to import ascii in R and store data matrix!"<<endl;
		return -1;
	}
	
	return 0;

}//close ImportDataInR()


int DataReader::ReadAscii(std::vector< std::vector<double> >& data,long int& ndim_read,std::string filename,std::string delimiter)
{
	//Open stream
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
    cout<<"ERROR: Cannot open input file "<<filename.c_str()<<"!"<<endl;
		return -1;
  }

	//Read ascii data line by line
	cout<<"INFO: Reading input data..."<<endl;
	std::string line;
	long int nentries= 0;	
	long int ndim= LONG_MAX;
	
	while(std::getline(in,line)) {
		
		char first_char= *(line.c_str());
		
		if(first_char!='#' && first_char!='\n' && first_char!=' '){

			//Parse current data line
			std::istringstream ss(line);
			std::vector<double> fields;
			std::string token;
			double value= 0;
			while(getline(ss,token,delimiter[0]))
			{
     		if(token=="nan" || token=="inf" || token=="-inf" || token=="na"){		
					value= TMath::SignalingNaN();
				}
				else {
					value= atof(token.c_str());
				}
				fields.push_back(value);
			}//end loop fields in this line

			//Append data to vector
			data.push_back(fields);
			long int ndim_curr= static_cast<long int>(fields.size());
			if(ndim_curr<ndim) ndim= ndim_curr;
			nentries++;

		}//close if

		if (!in.good()) break;

	}//close while events

	ndim_read= ndim;

	cout<<"INFO: #"<<nentries<<" entries read (ndim="<<ndim_read<<") ..."<<endl;

	//Close file
	in.close();

	return 0;

}//close ReadAscii()


TMatrixD* DataReader::ReadAscii(std::string filename,std::string delimiter)
{
	//Read ascii data
	long int ndim= 0;
	std::vector< std::vector<double> > data;
	if(ReadAscii(data,ndim,filename,delimiter)<0){
		cerr<<"ERROR: Failed to read data from file "<<filename<<"!"<<endl;
		return nullptr;
	}

	//Fill data matrix
	TMatrixD* dataMatrix= new TMatrixD(data.size(),ndim);
	for(size_t i=0;i<data.size();i++){
		for(long int j=0;j<ndim;j++){
			(*dataMatrix)(i,j)= data[i][j];
		}//end loop dim
	}//end loop entries

	return dataMatrix;

}//close ReadAscii()


}//close namespace

