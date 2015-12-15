// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double* k);
double ynt(const double teta, double* const y0, const double* const k, double dx);
double calct(const double* const y0, double* k, double dx);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
	ofstream out2("periode");
 	 const int dim = 2;
	double dx = 0.001,x=0;
	const double L = 100;
	  double y0[dim] = {1.01, 1.0};
	double yn[dim], k[8], teta;

	for (int i = 1; i <= 50; i++) {
	y0[0] = i/10.0;
	y0[1] = 0.0;
	x = 0;
  	out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << endl;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx, k);
		// check if we have passed a period 
	        if ((y0[1] > 0) and (yn[1] < 0)) {
        	        // perform root finding
			break;
        	}
    		for(int i=0; i<dim; i++) y0[i] = yn[i];
		out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << endl;
	}
        //cout << y0[1] << "\t" << yn[1] << endl;	
	teta = calct(y0, k, dx);
	out2 << i/10.0 << "\t" << x + teta*dx - dx << endl;
	cout << teta << endl;	
}
	out.close();
	out2.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k)
{
	const int dim = 2;
	double k1[dim], k2[dim], k3[dim], k4[dim];

 	 for(int i=0;i<dim; i++) {
	 k1[i] = y0[i];
	 k[i] = k1[i];
	}
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
 	 f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  	for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	
	k[2] = k2[0];
	k[3] = k2[1];
	k[4] = k3[0];
	k[5] = k3[1];
	k[6] = k4[0];
	k[7] = k4[1];	
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
	double y = y0[0];	
	
  	y0[0] = y0[1];
	y0[1] = -y/sqrt(1 + y * y);
}

double ynt(const double teta, const double* const y0, const double* const k, double dx) {
	//calculate bs
	double b1 = teta - (3 * teta * teta)/2.0 + (2 * teta * teta * teta)/3.0;
	double b2 = (teta * teta) - (2 * teta * teta * teta)/3.0;
	double b3 = b2;
	double b4 = -(teta * teta)/2.0 + (2 * teta * teta * teta)/3.0;
	
	return y0[1] + dx * (b1 * k[1] + b2 * k[3] + b3 * k[5] + b4 * k[7]); 	
}

double calct(const double* const y0, double* k, double dx) {
	double teta[3] = {0.0, 0.0, 1.0};
	double error = 1, ytemp;
	int counter = 0;
	while ((error > 1e-5) and (counter < 100)) {
		teta[1] = (teta[0] + teta[2])/2.0;

		ytemp = ynt(teta[1], y0, k, dx);
		if (ytemp  > 0) {
			teta[0] = teta[1];
		}		
		else {
			teta[2] = teta[1];
		}
		error = abs(ytemp);
		counter++;
		//cout <<counter <<  endl;
	}
	return teta[1];
}
