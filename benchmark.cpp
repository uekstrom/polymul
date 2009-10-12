#include <iostream>
#include <cmath>
#include "polymul.h"
#include "sys/time.h"

using namespace std;


template<class numtype, int Nvar, int Ndeg>
void benchmark(void)
{
  timeval t0,t1;
  double dt;
  polynomial<numtype, Nvar, Ndeg> p1,p2,pt;
  polynomial<numtype,Nvar,2*Ndeg> p3, sum;
  for (int i=0;i<p1.size;i++)
    {
      p1[i] = cos(i);
      p2[i] = cos(3*i);
    }
  sum = 0;

  cout << "Benchmarking " << Nvar  << ", " << Ndeg << ": " << 
    p1.size << " variables." << endl;

  // Multiplication
  int nrep = 2e7/p1.size;
  gettimeofday(&t0,NULL);
  for (int i=0;i<nrep;i++)
    {
      p1[0] = cos(i);
      polymul(p3,p1,p2);
      //   for (int j=0;j<sum.size;j++)
      //sum[j] += p3[j];
    }
  gettimeofday(&t1,NULL);
  dt = t1.tv_sec - t0.tv_sec + (t1.tv_usec - t0.tv_usec)*1e-6;
  cout << "polymul " << nrep/dt << " muls/s "  << sum[0] << " " <<  dt<< endl;

// Taylor Multiplication
  nrep = 20e7/p1.size;

  gettimeofday(&t0,NULL);
  for (int i=0;i<nrep;i++)
    {
      p1[0] = cos(i);
      taylormul(pt,p1,p2);
    }
  gettimeofday(&t1,NULL);
  dt = t1.tv_sec - t0.tv_sec + (t1.tv_usec - t0.tv_usec)*1e-6;
  cout << "taylormul " << nrep/dt << " muls/s "  << p1[0] << " " <<  dt<< endl;

  // Taylor in place Multiplication
  nrep = 20e6/p1.size;

  gettimeofday(&t0,NULL);
  for (int i=0;i<nrep;i++)
    {
      p1[0] = cos(i);
      taylormul(p1,p2);
    }
  gettimeofday(&t1,NULL);
  dt = t1.tv_sec - t0.tv_sec + (t1.tv_usec - t0.tv_usec)*1e-6;
  cout << "taylormul inplace " << nrep/dt << " muls/s "  << p1[0] << " " <<  dt<< endl;

  // Linear transformation
  nrep = 1e7/p1.size;
  numtype T[Nvar*Nvar];
  for (int i=0;i<Nvar*Nvar;i++)
    T[i] = sin(i);

  gettimeofday(&t0,NULL);
  for (int i=0;i<nrep;i++)
    {
      p1[0] = cos(i);
      T[0] = sin(i);
      polytrans(p1,p2,T);
    }
  gettimeofday(&t1,NULL);
  dt = t1.tv_sec - t0.tv_sec + (t1.tv_usec - t0.tv_usec)*1e-6;
  cout << "polytrans " << nrep/dt << " trans/s "  << p1[0] << " " <<  dt<< endl;
}

int main(void)
{
  benchmark<double,3,5>();
}
