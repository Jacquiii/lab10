#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u1, double* const u0, const double dx,const double dt, const double xmin,
                const int N);
void step(double* u2,double* u1,double* u0,
          const double dt, const double dx, const int N, double* h);

//---------------------------------------
int main(){

  const double tEnd = 0.2 ;


  const int N  = 1024;
  const double xmin = 0;
  const double xmax = 1;
  const double dx = (xmax-xmin)/N ;
  double dt = 0.001;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* u2 = new double[N];
  double* h;
  stringstream strm;

  initialize(u1,u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N);



  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      step(u2,u1,u0,dt,dx,N,h);
      
      h=u0;
      u0=u1;
      u1=u2;
      u2=h;
      
      t +=dt;
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
  }

  // Analytische Loesung
  initialize(u1,u0,dx,dt, xmin,N);
  for( int i=1; i<N-1; i++){
      double x = xmin + i*dx;
      double xi = x + u1[i]*tEnd;
      
      cout << xi << "\t" << u1[i] << endl; 
      
   }
  //writeToFile(u2, "Analytisch", dx, xmin,N);
  
  //cout << "t = " << t << endl;
  delete[] u0;
  delete[] u1;
  delete[] u2;
  return 0;
}
//-----------------------------------------------
void step(double* u2,double* u1,double* u0,
          const double dt, const double dx, const int N, double* h)
{
   u2[0] = u0[0] -dt/dx * u1[0] * (u1[1]-u1[N-1]); //Anfangsrandwert
   
   for( int i=1; i<N-1; i++){
      u2[i] = u0[i] - dt/dx * u1[i] * (u1[i+1]-u1[i-1]);
   }
   
   u2[N-1] = u0[N-1] - dt/dx * u1[N-1] * (u1[0]-u1[N-2]); //Endrandwert

   
}
//-----------------------------------------------
void initialize(double* const u1, double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{
   double u,ux, uxx;
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     u = sin(2*M_PI*x);
     ux = 2*M_PI*cos(2*M_PI*x);
     uxx = -4*pow(M_PI,2)*sin(2*M_PI*x);
     
     u1[i] = u;
     u0[i] = u -dt*ux + 1/2 * dt *uxx;
     
     
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << "\t" << endl;
   }
   out.close();
}
