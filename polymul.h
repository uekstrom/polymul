/*
Copyright (c) 2009 Ulf Ekstrom <uekstrom@gmail.com>

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef POLYMUL_H
#define POLYMUL_H
#include <cassert>

/*
  Multiply polynomials p1 and p2, both polynomials of nvar variables,
  of degree n1 and n2 respectively.  The result is a polynomial of
  degree n1+n2, which is _added_ to dst.

  A degree N polynomial of K variables has
  /N + K\
  \  N  /   terms.

  A polynomial P[K,N] is here organized as a sequence of monomials in
  K variables:

  M[K,0], M[K,1] .. M[K,N]  (for example 1, x y, x^2 xy y^2, ... )

  However, for each monomial all terms are of the same degree and the
  exponents of one variable can be inferred from the other
  exponents. We therefore have

  M[K,N] ~ P[K-1,N] (example: x^2 xy y^2 ~ 1 y y^2).

  To multiply monomials we can then use the ordinary polymul()
  function with lower nvar value.
*/

// Some compilers implement the C99 restrict also in C++.
#ifdef __GNUC__
#define POLYMUL_RESTRICT __restrict__
#else
#define POLYMUL_RESTRICT 
#endif


namespace polymul_internal
{

// Template to evaluate binomial coefficients at compile time.
template <int n, int k>
struct binomial
{
    enum { value = binomial<n-1,k-1>::value + binomial<n-1,k>::value };
};

template <int n>
struct binomial<n,n>
{
    enum { value = 1 };
};

template <int n>
struct binomial<n,0>
{
    enum { value = 1 };
};

template <int k>
struct binomial<0,k>
{
    enum { value = 1 };
};


// Recursive template classes for multiplication.

template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
class polynomial_multiplier
{
 public:
  // _add_ the product between p1 and p2 to dst.
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  { 
    polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2-1>::mul(dst,p1,p2);
    polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2>
      ::mul_monomial(dst,p1,p2+binomial<Nvar+Ndeg2-1,Ndeg2-1>::value);
  }
  // m2 is a monomial in Nvar variables and of order Ndeg2.
  // _add_ the product of p1 and m2 to dst.
  static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  { 
    polynomial_multiplier<numtype,Nvar-1,Ndeg1,Ndeg2>
      ::mul(dst+binomial<Nvar+Ndeg1+Ndeg2-1,Ndeg1+Ndeg2-1>::value,
	    p1 +binomial<Nvar+Ndeg1-1,Ndeg1-1>::value,
	    m2);
    polynomial_multiplier<numtype,Nvar,Ndeg1-1,Ndeg2>
      ::mul_monomial(dst,p1,m2);
  }


  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  { 
    polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2-1>
      ::mul_set(dst,p1,p2);
    // now we just add, because lower part is already set.
    polynomial_multiplier<numtype,Nvar,Ndeg1-1,Ndeg2>
      ::mul_monomial(dst,p1,p2+binomial<Nvar+Ndeg2-1,Ndeg2-1>::value);
    // Set final two highest monomials
    polynomial_multiplier<numtype,Nvar-1,Ndeg1,Ndeg2>
      ::mul_set(dst+ binomial<Nvar+Ndeg1+Ndeg2-1,Ndeg1+Ndeg2-1>::value,
		p1 + binomial<Nvar+Ndeg1-1,Ndeg1-1>::value,
		p2 + binomial<Nvar+Ndeg2-1,Ndeg2-1>::value);
  }

  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  { 
    // Unclear what this would mean, since the lower parts of dst are not touched.
    assert(0 && " I am not supposed to be here..");
  }

};

template<class numtype, int Ndeg1, int Ndeg2>
  class polynomial_multiplier<numtype, 1, Ndeg1, Ndeg2>
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<=Ndeg1;i++)
      for (int j=0;j<=Ndeg2;j++)
	dst[i+j] += p1[i]*p2[j];
  }
  static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=0;i<=Ndeg1;i++)
      dst[i+Ndeg2] += m2[0]*p1[i];
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<=Ndeg1;i++)
      dst[i] = p1[i]*p2[0]; 
    // Now set all where i = Ndeg1 and j > 0
    for (int j=1;j<=Ndeg2;j++)
      dst[Ndeg1+j] = p1[Ndeg1]*p2[j];
    // Add all inbetween
    for (int i=0;i<Ndeg1;i++)
      for (int j=1;j<=Ndeg2;j++)
	dst[i+j] += p1[i]*p2[j];
  }
  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=0;i<Ndeg1;i++)
      dst[i+Ndeg2] = m2[0]*p1[i];
  }
};

template<class numtype, int Nvar, int Ndeg2>
  class polynomial_multiplier<numtype, Nvar, 0, Ndeg2>
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<Nvar+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*p2[i];
  }
  static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=binomial<Nvar+Ndeg2-1,Ndeg2-1>::value;i<binomial<Nvar+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*m2[i-binomial<Nvar+Ndeg2-1,Ndeg2-1>::value];
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<Nvar+Ndeg2,Ndeg2>::value;i++)
      dst[i] = p1[0]*p2[i];
  }
  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=binomial<Nvar+Ndeg2-1,Ndeg2-1>::value;i<binomial<Nvar+Ndeg2,Ndeg2>::value;i++)
      dst[i] = p1[0]*m2[i-binomial<Nvar+Ndeg2-1,Ndeg2-1>::value];
  }
};

template<class numtype, int Ndeg2>
  class polynomial_multiplier<numtype, 1, 0, Ndeg2>
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<1+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*p2[i];
  }
  static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=binomial<1+Ndeg2-1,Ndeg2-1>::value;i<binomial<1+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*m2[i-binomial<1+Ndeg2-1,Ndeg2-1>::value];
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<1+Ndeg2,Ndeg2>::value;i++)
      dst[i] = p1[0]*p2[i];
  }
  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=binomial<1+Ndeg2-1,Ndeg2-1>::value;i<binomial<1+Ndeg2,Ndeg2>::value;i++)
      dst[i] = p1[0]*m2[i-binomial<1+Ndeg2-1,Ndeg2-1>::value];
  }
};

template<class numtype, int Nvar, int Ndeg1>
  class polynomial_multiplier<numtype, Nvar, Ndeg1, 0>
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<Nvar+Ndeg1,Ndeg1>::value;i++)
      dst[i] += p1[i]*p2[0];
  }
  static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=0;i<binomial<Nvar+Ndeg1,Ndeg1>::value;i++)
      dst[i] += p1[i]*m2[0];
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<Nvar+Ndeg1,Ndeg1>::value;i++)
      dst[i] = p1[i]*p2[0];
  }
  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=0;i<binomial<Nvar+Ndeg1,Ndeg1>::value;i++)
      dst[i] = p1[i]*m2[0];
  }
};

template<class numtype, int Nvar>
  class polynomial_multiplier<numtype, Nvar, 0, 0>
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    dst[0] += p1[0]*p2[0];
  }
  static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    dst[0] += p1[0]*m2[0];
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    dst[0] = p1[0]*p2[0];
  }
  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {
    dst[0] = p1[0]*m2[0];
  }
};


// Like polymul but truncates the result to Ndeg1 order.
// Ndeg2 _must_ be less or equal to Ndeg1, otherwise the recursion
// will not terminate.
template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
class taylor_multiplier
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    polynomial_multiplier<numtype,Nvar,Ndeg1-Ndeg2,Ndeg2>
      ::mul_monomial(dst,p1,p2+binomial<Nvar+Ndeg2-1,Ndeg2-1>::value);   
    taylor_multiplier<numtype,Nvar,Ndeg1,Ndeg2-1>::mul(dst,p1,p2); 
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    taylor_multiplier<numtype,Nvar,Ndeg1,Ndeg2-1>::mul_set(dst,p1,p2); 
    polynomial_multiplier<numtype,Nvar,Ndeg1-Ndeg2,Ndeg2>
      ::mul_monomial(dst,p1,p2+binomial<Nvar+Ndeg2-1,Ndeg2-1>::value);   
  }
};

template<class numtype, int Nvar, int Ndeg1>
  class taylor_multiplier<numtype, Nvar,Ndeg1,0>
{
 public:
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<Ndeg1+Nvar,Ndeg1>::value;i++)
      dst[i] += p1[i]*p2[0];
  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<binomial<Ndeg1+Nvar,Ndeg1>::value;i++)
      dst[i] = p1[i]*p2[0];
  }
};


template<class numtype, int Nvar, int Ndeg, int i2> // (2)
class taylor_inplace_multiplier
{
 public:
  static void mul(numtype p1[], const numtype p2[])
  {
    // M1(Ndeg-i2)*M2(i2) -> M1(Ndeg)
    polynomial_multiplier<numtype,Nvar-1,Ndeg-i2,i2>
      ::mul(p1+binomial<Nvar+Ndeg-1,Ndeg-1>::value,
	    p1+binomial<Nvar+Ndeg-i2-1,Ndeg-i2-1>::value,
	    p2+binomial<Nvar+i2-1,i2-1>::value);
    taylor_inplace_multiplier<numtype,Nvar,Ndeg,i2+1> // Back to (2), or to (3) when i2+1 == Ndeg
      ::mul(p1,p2);
  }
};

template<class numtype, int Nvar, int Ndeg>
class taylor_inplace_multiplier<numtype, Nvar, Ndeg, Ndeg> // (3) final contribution to M1(Ndeg)
{
 public:
  static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
  {
    for (int i=binomial<Nvar+Ndeg-1,Ndeg-1>::value;
	 i<binomial<Nvar+Ndeg,Ndeg>::value;i++)
      p1[i] += p2[i]*p1[0];
    taylor_inplace_multiplier<numtype,Nvar,Ndeg-1,0>::mul(p1,p2); // Do lower degree terms, or go to (4)
  }
};

template<class numtype, int Nvar, int Ndeg>
class taylor_inplace_multiplier<numtype, Nvar, Ndeg, 0> // (1) comes in here, sets M1(Ndeg)
{
 public:
  static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
  {
    for (int i=binomial<Nvar+Ndeg-1,Ndeg-1>::value;i<binomial<Nvar+Ndeg,Ndeg>::value;i++)
      p1[i] *= p2[0];
    taylor_inplace_multiplier<numtype,Nvar,Ndeg,1>::mul(p1,p2); // continue at (2)
  }
};

template<class numtype, int Nvar>
class taylor_inplace_multiplier<numtype, Nvar, 0, 0> // (4), last coefficient.
{
 public:
  static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
  {
    p1[0] *= p2[0];
  }
};

template<class numtype>
class taylor_inplace_multiplier<numtype, 1, 0, 0> // (4), last coefficient.
{
 public:
  static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
  {
    p1[0] *= p2[0];
  }
};

template<class numtype, int Ndeg>
class taylor_inplace_multiplier<numtype, 1, Ndeg, 0> //Above code is only for Ndeg>1
{
 public:
  static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
  {
    p1[Ndeg] *= p2[0];
    for (int i=0;i<Ndeg;i++)
      p1[Ndeg] += p1[i]*p2[Ndeg-i];
    taylor_inplace_multiplier<numtype, 1, Ndeg-1, 0>::mul(p1,p2);
  }
};

 template<class numtype, class vartype, int Nvar, int Ndeg>
class polynomial_evaluator
{
 public:
  static vartype eval(const numtype p[], const vartype x[])
  {
    return polynomial_evaluator<numtype,vartype,Nvar,Ndeg>
      ::eval_monomial(p+binomial<Nvar+Ndeg-1,Ndeg-1>::value,x) +
      polynomial_evaluator<numtype,vartype,Nvar,Ndeg-1>
      ::eval(p,x);
  }

  // Evaluate monomial in Nvar variables.
  // M(N,K) = x[0]*M(N-1,K) + M(N,K-1)
  static vartype eval_monomial(const numtype p[], const vartype x[])
  {
    return x[0]*polynomial_evaluator<numtype,vartype,Nvar,Ndeg-1>
      ::eval_monomial(p,x) +
      polynomial_evaluator<numtype,vartype,Nvar-1,Ndeg>
      ::eval_monomial(p+binomial<Nvar+Ndeg-2,Ndeg-1>::value,
		      x+1);
  }
};


/*
  1 x y x^2 xy y^2 =
  
  1 [x] [y] x[x y] y^2

 */
template<class numtype, class vartype, int Ndeg>
class polynomial_evaluator<numtype, vartype, 1, Ndeg>
{
 public:
  // a + bx + cx^2 = a + x(b + x(c))
  static vartype eval(const numtype p[], const vartype x[])
  {
    // Horner scheme:
    vartype sum = p[Ndeg];
    for (int i=Ndeg-1;i>=0;i--)
      sum = sum*x[0] + p[i];
    return sum;
  }
  static vartype eval_monomial(const numtype p[], const vartype x[])
  {
    vartype xn = 1;
    for (int i=0;i<Ndeg;i++)
      xn *= x[0];
    return xn*p[0];
  }
};

template<class numtype, class vartype, int Nvar>
class polynomial_evaluator<numtype, vartype, Nvar, 0>
{
 public:
  static vartype eval(const numtype p[], const vartype x[])
  {
    return p[0];
  }
  static vartype eval_monomial(const numtype p[], const vartype x[])
  {
    return p[0];
  }
};

}
// End of namespace polymul_internal


template<class numtype, int Nvar, int Ndeg>
  class polynomial
{
 public:
  polynomial(void) {}
  polynomial(const numtype &c0) 
    { 
      c[0] = c0;
      for (int i=1;i<this->size();i++)
	c[i] = 0;
    }
  numtype operator[](int i) const
  {
    assert(i>=0);
    assert(i<this->size());
    return c[i];
  }
  numtype &operator[](int i)
  {
    assert(i>=0);
    assert(i<this->size());
    return c[i];
  }
  static int size(void) { return polymul_internal::binomial<Nvar+Ndeg,Ndeg>::value; }
  void zero(void)
  {
    for (int i=0;i<size();i++)
      c[i] = 0;
  }
  // This is a _very slow_ function to get the exponents
  // of a particular term. 
  static void exponents(int term, int exponents[Nvar])
  {
    assert(term >= 0);
    if (Nvar == 1)
      {
	exponents[0] = term;
	return;
      }
    for (int i=0;i<Nvar;i++)
      exponents[i] = 0;
    if (term >= polymul_internal::binomial<Nvar+Ndeg,Ndeg>::value)
      {
	assert(0 && "term < binomial<Nvar+Ndeg,Ndeg>::value");
      }
    for (int i=0;i<term;i++)
      polynomial<numtype,Nvar,Ndeg>::next_exponents(Nvar,exponents);
  }
  // Return the index of the term with certain exponents
  static int term_index(const int exponents[Nvar])
  {
    int N = 0;
    for (int i=0;i<Nvar;i++)
      N += exponents[i];
    int i = 0, idx = 0;
    N--;
    while (N >= 0)
      {
	idx += polynomial<numtype,Nvar,Ndeg>::polylen(Nvar-i,N);
	N -= exponents[i];
	i++;
      }
    return idx;
  }
  // Evaluate the polynomial at x.
  template<class vartype>
  numtype eval(const vartype x[Nvar]) const
  {
    return polymul_internal::polynomial_evaluator<numtype,vartype,Nvar,Ndeg>::
      eval(c,x);
  }

  numtype c[polymul_internal::binomial<Nvar+Ndeg,Ndeg>::value];

 protected:
  static void next_exponents(int nvar, int m[Nvar])
  {
    int k = 0;
    for (int i=0;i<nvar-1;i++)
      k += m[i];
    if (k == 0)
      {
	m[0] = m[nvar-1] + 1;
	m[nvar-1] = 0;
	return;
      }
    if (m[nvar-2] > 0)
      {
	m[nvar-1]++;
	m[nvar-2]--;
      }
    else
      {
	next_exponents(nvar-1,m);
	for (int i=nvar-2;i>=0;i--)
	  {
	    if (m[i] > 0)
	      {
		m[i] += m[nvar-1];
		break;
	      }
	  }
	m[nvar-1] = 0;
      }
  }
  // = binomial(nvar+ndeg,ndeg), but this one
  // can be evaluated at run time.
  static int polylen(int nvar, int ndeg)
  {
    int len = 1;
    for (int k=1;k<=nvar;k++)
      {
	len *= ndeg + k;
	len /= k;
      }
    return len;
  }
};


// User interface

template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
  void polymul(polynomial<numtype, Nvar,Ndeg1+Ndeg2> & POLYMUL_RESTRICT dst,
	       const polynomial<numtype, Nvar,Ndeg1> &p1,
	       const polynomial<numtype, Nvar,Ndeg2> &p2)
{
  polymul_internal::polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2>
    ::mul_set(dst.c,p1.c,p2.c);
}

template<class numtype, int Nvar, int Ndeg>
  void taylormul(polynomial<numtype, Nvar,Ndeg> & POLYMUL_RESTRICT dst,
		 const polynomial<numtype, Nvar,Ndeg> &p1,
		 const polynomial<numtype, Nvar,Ndeg> &p2)
{
  polymul_internal::taylor_multiplier<numtype,Nvar,Ndeg,Ndeg>
    ::mul_set(dst.c,p1.c,p2.c);
}

template<class numtype, int Nvar, int Ndeg>
void taylormul(polynomial<numtype, Nvar,Ndeg> & POLYMUL_RESTRICT p1,
	       const polynomial<numtype, Nvar,Ndeg> &p2)
{
  polymul_internal::taylor_inplace_multiplier<numtype,Nvar,Ndeg,0>
    ::mul(p1.c,p2.c);
}

#endif
