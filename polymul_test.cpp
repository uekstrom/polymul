#include <iostream>
#include <valarray>
using namespace std;
#include "polymul.h"



template<class num, int Nvar, int Ndeg>
void printpoly(ostream &dst, const polynomial<num,Nvar,Ndeg> &p)
{
  int exps[Nvar] = {0};
  for (int i=0;i<p.size();i++)
    {
      p.exponents(i,exps);
      for (int j=0;j<Nvar;j++)
	cout << exps[j] << " ";
      cout << p[i] << endl;
      assert(i == p.term_index(exps));
    }
}

int main(void)
{
  int pp_good[] = 
    {0, -7, 0, 7, 20, 40, 82, 119, 177, 227, 21, 
     80, 134, 185, 367, 365, 153, 400, 534, 326, 
     44, 125, 184, 327, 559, 561, 302, 782, 976, 
     595, 238, 591, 1066, 849, 496};

  int eval_1d_x[] = {11,-13,17};
  int eval_1d_p[] = {3,5,-7};
  int eval_1d_px[] = {-789, -1245, -1935};

  cout << "Starting tests.." << endl;
  

  // Evaluating
  polynomial<int,1,2> peval_1d;
  for (int i=0;i<peval_1d.size();i++)
    peval_1d[i] = eval_1d_p[i];
  for (int i=0;i<3;i++)
    if (peval_1d.eval(&eval_1d_x[i]) != eval_1d_px[i])
      {
	cout << "WARNING (polyval 1d), expected " << eval_1d_px[i] <<
	  " got " << peval_1d.eval(&eval_1d_x[i]) << endl;
      }
  
  // Multiplying
  polynomial<int,3,2> p1, p2, pt;
  polynomial<int,3,4> pp(1);
  for (int i=0;i<p1.size();i++)
    p1[i] = 7+i;
  for (int i=0;i<p2.size();i++)
    p2[i] = i*i/2 - i;
  
  polymul(pp,p1,p2);
  for (int i=0;i<pp.size();i++)
    {
      if (pp_good[i] != pp[i])
	{
	  cout << "WARNING (known good vs polymul): " <<  i << " "<< pp_good[i]<< " " << pp[i] <<endl;
	}
    }
  // Check that multiply then eval gives the same as eval then multiply
  int x1[3] = {12,-1,7};
  int x2[3] = {-4,3,1};
  if (p1.eval(x1)*p2.eval(x1) != pp.eval(x1) || 
      p1.eval(x2)*p2.eval(x2) != pp.eval(x2))
    {
      cout << "WARNING: multidimensional eval failed\n";
    }

  //start of polymul result should match taylormul
  taylormul(pt,p2,p1);
  for (int i=0;i<pt.size();i++)
    if (pt[i] != pp[i])
      {
	cout << "WARNING (taylor vs poly): " <<  i << " "<< pp[i]<< " " << pt[i] <<endl;
      }

  taylormul(p1,p2);
  for (int i=0;i<pt.size();i++)
    {
      if (p1[i] != pt[i])
	{
	  cout <<  i << " "<< p1[i]<< " " << pt[i] << " WARNING (in-place)" << endl;
	}
    }
  cout << "If no warnings were printed above, then things are fine." << endl;
  cout << "Polynomial exponents and coefficients:\n";
  printpoly(cout,p1);
  cout << "End of tests." << endl;
  return 0;
}
