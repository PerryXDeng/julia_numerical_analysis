// Interpolation: Splines
//
// J. Faber, for MATH-411 (RIT)
// GPL, etc.
//
// call it with ./spline

//You need an input file spline_data.dat
//with the first line containing the number of points
//the second line the spline method 
//the third lines is the values of the endpoint slopes, which will be ingored if nsp=1,4, or 5; the line
//is still necessary!!!
// and then the x y values, like this.  First and last are endpoints.

// Spline methods are
// 1: Natural splines: y_1'' = y_n'' = 0
// 2: Curvature adjusted: PASS IN VALUES of y_1'' and y_n''
// 3: Clamped: PASS IN VALUES of y_1' and y_n'
// 4: parabolically terminated: y_1''' = y_n''' = 0
// 5: not-a-knot: y_1''' = y_2'''; y_n''' = y_{n-1}'''

/*
6
2 
3 4 
0.5 0
1.0 1.0
2.0 0.0
3.0 -1.0
4.0 4.0
5.5 3.0
*/
// Interpolated values are output to spline_out.dat

// These headers are so standard one tends not to even think about them
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "stdio.h"
#include "stdlib.h"

//THIS IS A GSL_BASED CODE; compile with at least -lgsl, plus links to libraries and include files as necessary
#include <gsl/gsl_linalg.h>

using namespace std;

int main(int argc, char* argv[]) {

  //We'll read in the data from a file
  ifstream infile;
  infile.open("spline_data.dat");
  ofstream outfile;
  outfile.open("spline_out.dat");

  //First line of the file should be thenumber of points
  int n;
  infile>>n;


  //Now we read in the spline method and two paameters
  int nsp;
  double c1,cn;

  infile>>nsp;
  infile>>c1>>cn;

  //then the data values (x_i,y_i);
  double* xpt=new double[n];
  double* ypt=new double[n];
  int i,j;

  for (i=0; i<n; i++){ infile>>xpt[i]>>ypt[i];
    cout<<i<<" "<<xpt[i]<<" "<<ypt[i]<<endl;}

  double* matrix = new double[n*n];
  double* RHS = new double[n];
  double* solution = new double[n];
	
  for (i=0; i<n; i++) {
    RHS[i]=0.;
    solution[i]=0.;
    for (j=0; j<n; j++) {
      matrix[i*n+j]=0.;
    }
  }
	
  for (i=1; i<n-1; i++) {
    
    RHS[i] = 6.0*((ypt[i+1]-ypt[i])/(xpt[i+1]-xpt[i])-
		  (ypt[i]-ypt[i-1])/(xpt[i]-xpt[i-1]));

    for (j=0; j<n; j++) {
      int index = i*n+j;
      
      //Three components for splines
      if(j==i-1) matrix[index]+=xpt[i]-xpt[i-1];
      if(j==i) matrix[index]+=2.0*(xpt[i+1]-xpt[i-1]);
      if(j==i+1) matrix[index]+=xpt[i+1]-xpt[i];
    }
  }
  
  //Natural spline Boundary conditions
  if(nsp==1) {
    matrix[0] = 1.0;
    RHS[0]=0.0;
    matrix[n*n-1] = 1.0;
    RHS[n-1] = 0.0;

  } else if (nsp==2) {
    //Curvature-adjusted
    matrix[0] = 1.0;
    RHS[0] = c1;
    matrix[n*n-1] = 1.0;
    RHS[n-1] = cn;

  } else if (nsp==3) {
    //Clamped
    matrix[0] = 2.0;
    matrix[1] = 1.0;
    RHS[0]=6.0*((ypt[1]-ypt[0])/(xpt[1]-xpt[0]) - c1);
    matrix[n*n-2] = 1.0;
    matrix[n*n-1] = 2.0;
    RHS[n-1] = 6.0*(cn - (ypt[n-1]-ypt[n-2])/(xpt[n-1]-xpt[n-2]) );

  }  else if (nsp==4) {
    //parabolically terminated
    matrix[0] = 1.0;
    matrix[1] = -1.0;
    RHS[0]=0.0;
    matrix[n*n-2] = -1.0;
    matrix[n*n-1] = 1.0;
    RHS[n-1] = 0.0;

  }   else if (nsp==5) {
    //not-a-knot
    matrix[0] = xpt[2]-xpt[1];
    matrix[1] = xpt[0]-xpt[2];
    matrix[2] = xpt[1]-xpt[0];
    RHS[0]=0.0;
    matrix[n*n-3] = xpt[n-1]-xpt[n-2];
    matrix[n*n-2] = xpt[n-3]-xpt[n-1];
    matrix[n*n-1] = xpt[n-2]-xpt[n-3];
    RHS[n-1] = 0.0;
  } 
    
  //This is code to call the gnu scientific library.  It solved Ax=b in Matrix form using
  //an LU decomposition with partial pivoting.
  //See www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html for the source

  gsl_matrix_view m = gsl_matrix_view_array (matrix,n,n);
  gsl_vector_view b = gsl_vector_view_array (RHS,n);
  gsl_vector *x = gsl_vector_alloc(n);
  int s;
  gsl_permutation * p = gsl_permutation_alloc(n);
  gsl_linalg_LU_decomp(&m.matrix,p,&s);
  gsl_linalg_LU_solve(&m.matrix,p,&b.vector,x);

  for (i=0; i<n; i++) {
    solution[i] = gsl_vector_get(x, i);
    cout<<i<<" "<<solution[i]<<endl;
  }

  gsl_permutation_free(p);
  gsl_vector_free(x);


  //And back to our regularly scheduled program.  We will output stuff!

  int ninterval=10;

  outfile<<xpt[0]<<" "<<ypt[0]<<endl;

  //Having output the beginning point, loop over points
  for (i=0; i<n-1; i++) {
    double dx = xpt[i+1]-xpt[i];
    
    for (j=1; j<=ninterval; j++) {

      //These are the A and B values for the spline
      double bb = 1.0*j/(1.0*ninterval);
      double aa=1.0-bb;
      double xval = xpt[i]+bb*dx;

      double cc = dx*dx*aa*(aa*aa-1.0)/6.0;
      double dd = dx*dx*bb*(bb*bb-1.0)/6.0;
      
      double yval = aa*ypt[i]+bb*ypt[i+1]+cc*solution[i]+dd*solution[i+1];
      
      outfile<<xval<<" "<<yval<<endl;
    }
  }

  outfile.close();
  
  return 0;
  
}

