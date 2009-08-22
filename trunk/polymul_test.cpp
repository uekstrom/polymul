#include <iostream>
#include "polymul.h"

using namespace std;

template<class num, int Nvar, int Ndeg>
void printpoly(ostream &dst, const polynomial<num,Nvar,Ndeg> &p)
{
  int exps[Nvar] = {0};
  for (int i=0;i<p.size;i++)
    {
      p.exponents(i,exps);
      for (int j=0;j<Nvar;j++)
	cout << exps[j] << " ";
      cout << " " << p[i] << endl;
      assert(i == p.term_index(exps));
    }
}

#ifdef __GNUC__NOT_ICC__
// Test polynomials over gcc vectors
typedef double v2d __attribute__ ((vector_size (16)));
void gcc_vectest(void)
{
  v2d v0 = {0,0}, c0 = {1,2};
  polynomial<v2d,2,2> p1,p2,pp;
  for (int i=0;i<p1.size;i++)
    {
      p1[i] = v0;
      p2[i] = v0;
    }
  p1[0] = c0;
  p2[1] = c0;
  taylormul(pp,p1,p2);
}
#endif

int fac(int n)
{
  assert(n>=0);
  assert(n<13);
  int f = 1;
  while(n>1)
    {
      f *= n;
      n--;
    }
  return f;
}

int main(void)
{
  int pp_good[] = 
    {0, -7, 0, 7, 20, 40, 82, 119, 177, 227, 21, 
     80, 134, 185, 367, 365, 153, 400, 534, 326, 
     44, 125, 184, 327, 559, 561, 302, 782, 976, 
     595, 238, 591, 1066, 849, 496};

  int eval_1d_x[] = {11,-13,17};
  double eval_1d_x_double[] = {11,-13,17};
  int eval_1d_p[] = {3,5,-7};
  int eval_1d_px[] = {-789, -1245, -1935};

  cout << "Starting tests.." << endl;
  

  // Evaluating
  polynomial<int,1,2> peval_1d;
  for (int i=0;i<peval_1d.size;i++)
    peval_1d[i] = eval_1d_p[i];
  for (int i=0;i<3;i++)
    if (peval_1d.eval(&eval_1d_x_double[i]) != eval_1d_px[i])
      {
	cout << "WARNING (polyval 1d), expected " << eval_1d_px[i] <<
	  " got " << peval_1d.eval(&eval_1d_x_double[i]) << endl;
      }
  
  // Multiplying
  polynomial<int,3,2> p1, p2, pt;
  polynomial<int,3,4> pp(1);
  for (int i=0;i<p1.size;i++)
    p1[i] = 7+i;
  for (int i=0;i<p2.size;i++)
    p2[i] = i*i/2 - i;
  
  polymul(pp,p1,p2);
  for (int i=0;i<pp.size;i++)
    {
      if (pp_good[i] != pp[i])
	{
	  cout << "WARNING (known good vs polymul): " <<  i << " "<< pp_good[i]<< " " << pp[i] <<endl;
	}
    }
  // Is term_prod correct?
  if (polymul_internal::term_prod<3,1,2>::prod != 5)
    cout << "WARNING: term_prod<3,1,2> incorrect = " << polymul_internal::term_prod<3,1,2>::prod << endl;
  if (polymul_internal::term_prod<3,3,3>::prod != 9)
    cout << "WARNING: term_prod<3,3,3> incorrect = " << polymul_internal::term_prod<3,3,3>::prod << endl;

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
  for (int i=0;i<pt.size;i++)
    if (pt[i] != pp[i])
      {
	cout << "WARNING (taylor vs poly): " <<  i << " "<< pp[i]<< " " << pt[i] <<endl;
      }

  taylormul(p1,p2);
  for (int i=0;i<pt.size;i++)
    {
      if (p1[i] != pt[i])
	{
	  cout <<  i << " "<< p1[i]<< " " << pt[i] << " WARNING (in-place)" << endl;
	}
    }
  // taylormul with lower order p2
  polynomial<int,4,5> pleft,ppadded(0),pcheck1,pcheck2;
  polynomial<int,4,2> plower;
  for (int i=0;i<plower.size;i++)
    {
      plower[i] = i+1;
      ppadded[i] = plower[i];
    }
  for (int i=0;i<pleft.size;i++)
    pleft[i] = 1-i+((i*i) % 9);
  pcheck1 = pleft;
  taylormul(pcheck1,ppadded);
  pcheck2 = pleft;
  taylormul(pcheck2,plower);
  for (int i=0;i<pcheck1.size;i++)
    {
      if (pcheck1[i] != pcheck2[i])
	cout << "taylormul lower order error at " << i << ", expected " 
	     << pcheck1[i] << " got " << pcheck2[i] << endl;
    }
  
  // Test contraction
  // dot with 1 2 3 ..
  polynomial<int,3,4> pd;
  polynomial<int,3,2> pd1;
  for (int i=0;i<pd.size;i++)
    pd[i] = i+1;
  double dot1 = 0, dot2 = 0;
  // Reuse p1 and p2 from above
  polymul(pp,p1,p2);
  for (int i=0;i<pp.size;i++)
    dot1 += pp[i]*pd[i];
  polycontract(pd1,p2,pd);
  for (int i=0;i<p1.size;i++)
    dot2 += p1[i]*pd1[i];
  if (dot1 != dot2)
    {
      cout << "WARNING, in contraction: dot1 and dot2: " << dot1 << " " << dot2 << endl;
    }

  if (polymul_internal::term_deg<3,4>::deg != 2 or
      polymul_internal::term_deg<2,6>::deg != 3)
    {
      cerr << "WARNING, error in term_deg<>" << endl;
    }
  // Sparse (single term) multiply
  p1 = 0;
  p2 = 0;
  p1[2] = 3;
  p2[4] = 7;
  polymul(pp,p1,p2);
  polynomial<int,3,4> ppz;
  ppz = 0;
  polymul_term<int,3,2,4>(ppz,p1,7);
  for (int i=0;i<ppz.size;i++)
    if (pp[i] != ppz[i])
      cout << "WARNING: single term multiply failed at " << i << endl;


#if 0
  //Eval terms
  double ex[] = {2,3,5,1,1};
  polynomial<double,5,5> pterms;
  pterms.zero();
  polyterms(pterms,ex);
  cout << "Terms for x=2, y=3, z=5:\n";
  printpoly(cout,pterms);
#endif
#ifdef BENCHMARK
  // Eval bench
  double sum = 0;
  for (int i=0;i<1e7;i++)
    {
      ex[1] += i & 2;
      sum += pterms.eval_new(ex);
      sum /= (sum - 1);
    }
  cout << sum << endl;
#endif
#ifdef __GNUC__NOT_ICC__
  gcc_vectest();
#endif  
  cout << "If no warnings were printed above, then things are fine." << endl;
  cout << "Polynomial exponents and coefficients:\n";
  printpoly(cout,p1);
  cout << "End of tests." << endl;
  return 0;
}
