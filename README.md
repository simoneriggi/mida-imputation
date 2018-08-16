<p align="left">
  <img src="share/logo.png" alt="Sample outputs"/>
</p>

# mida-imputation
A C++ software tool for imputating multivariate missing data with multiple methods (mean, multiple imputation, least-wise deletion, multivariate normal and skew-normal clustering). It is distributed for research use only under the GNU General Public License v3.0.

## **Credits**
If you use this software for your research, please acknowledge it in your papers by citing the following references:

* `S. Riggi et al., "Combined spectrum-Xmax fit with a likelihood approach", Nuclear Instruments and Methods in Physics Research A 780 (2015) 81–90`

or consider including me (`S. Riggi, INAF - Osservatorio Astrofisico di Catania, Via S. Sofia 78, I-95123, Catania, Italy`)
as a co-author on your publications.

## **Status**
Software is currently been updated.

## **Installation**  

### **Prerequisites**
Install the project mandatory dependencies:  
* ROOT [https://root.cern.ch/]
* R [https://www.r-project.org/], install also these additional packages: RInside, Rcpp, Matrix, Amelia, flexclust

Make sure you have set the following environment variables to the external library installation dirs 
* ROOTSYS: set to ROOT installation prefix

NB: Modify Makefile CPPFLAGS and LDFLAGS in case the dependency tools cannot be found.

### **Build**
To build the project:

* Clone this repository into your local $SOURCE_DIR    
  ```git clone https://github.com/simoneriggi/mida-imputation.git $SOURCE_DIR```    
* In the project directory type:    
  ```make```  

Binaries will be placed in the bin/ directory and libraries in the lib/ directory.

### **Usage**
* ```MDImputation [--input=[path-to-inputfile]] [--config=[path-to-configfile]]```    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```--input=[path-to-inputfile] -  Input data file (.dat) with missing data to be imputed```   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```--input=[path-to-inputfile] -  Input data file (.dat) with missing data to be imputed```   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```--method=[imputation-method] -  Imputation method to be used (1=MEAN, 2=LISTWISE DELETION, 3=MultipleImputation, 4=MN clustering, 5=MSN clustering```   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```--config=[path-to-configfile] - Configuration file name with options```    
