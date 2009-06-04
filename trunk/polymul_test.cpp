#include <iostream>
#include "polymul.h"

using namespace std;

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
  cout << "Starting tests.." << endl;

  polynomial<int,3,2> p1, p2, pt;
  polynomial<int,3,4> pp(0),ppset;
  for (int i=0;i<p1.size();i++)
    p1[i] = 7+i;
  for (int i=0;i<p2.size();i++)
    p2[i] = i*i/2 - i;

  pp.zero();
  polymul_add(pp,p1,p2);
  polymul(ppset,p1,p2);
  for (int i=0;i<pp.size();i++)
    {
      if (pp_good[i] != pp[i])
	{
	  cout << "WARNING (known good vs polymul_add): " <<  i << " "<< pp_good[i]<< " " << pp[i] <<endl;
	}
      if (ppset[i] != pp[i])
	{
	  cout << "WARNING (set vs add): " <<  i << " "<< pp[i]<< " " << ppset[i] <<endl;
	}
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
