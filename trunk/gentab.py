# Program to generate specific template instances for most of the operations
# that the polymul library can do. This may give a significal speedup, but can
# also affec the size of the generated code. It all depends on the compiler,
# really. By Ulf Ekstrom <uekstrom@gmail.com> 2009

def exponents(nvar, ndeg):
    """Return a list of all exponents of an (nvar,ndeg) polynomial, in graded lexicographical order"""
    if ndeg == 0:
        return [ [0]*nvar ]
    elif nvar == 0:
        return [[]]
    else:
        l = []
        for n in range(ndeg+1):
            lower = exponents(nvar-1,n)
            for k in lower:
                l.append([n-sum(k)]+k)
        return l


def monomial_exponents(nvar, ndeg):
    return [e for e in exponents(nvar,ndeg) if sum(e) == ndeg]
        
def expadd(t1,t2):
    return [a[0]+a[1] for a in zip(t1,t2)]

def taylormul_instance(nvar, ndeg1, ndeg2, maxterms = -1):
    terms1 = exponents(nvar,ndeg1)
    terms2 = exponents(nvar,ndeg2)
    terms_res = terms1
    prods = []
    for (i1,t1) in enumerate(terms1):
        for (i2,t2) in enumerate(terms2):
            prodt = expadd(t1,t2)
            if sum(prodt) <= ndeg1:
                prods.append((terms_res.index(prodt),i1,i2))
    if maxterms > -1 and len(prods) > maxterms:
        return False
    prods.sort()
    print """template<class numtype>
struct taylor_multiplier<numtype, %i, %i, %i>
{
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {""" % (nvar,ndeg1,ndeg2)
    for k in range(len(terms_res)):
#    for k in range(len(terms_res)-1,-1,-1):
        pk = ["p1[%i]*p2[%i]" % (p[1],p[2]) for p in prods if p[0] == k]
        print "    dst[%i] +=" % k," + ".join(pk),";"
    print """  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {"""
    for k in range(len(terms_res)):
#    for k in range(len(terms_res)-1,-1,-1):
        pk = ["p1[%i]*p2[%i]" % (p[1],p[2]) for p in prods if p[0] == k]
        print "    dst[%i] =" % k," + ".join(pk),";"
    print """  }
};"""
    return True
            

def taylormul_inplace_instance(nvar, ndeg1, ndeg2, maxterms = -1):
    terms1 = exponents(nvar,ndeg1)
    terms2 = exponents(nvar,ndeg2)
    terms_res = terms1
    prods = []
    for (i1,t1) in enumerate(terms1):
        for (i2,t2) in enumerate(terms2):
            prodt = expadd(t1,t2)
            if sum(prodt) <= ndeg1:
                prods.append((terms_res.index(prodt),i1,i2))
    if maxterms > -1 and len(prods) > maxterms:
        return False
    prods.sort()
    print """template<class numtype>
struct taylor_inplace_multiplier<numtype, %i, %i, %i, 0>
{
  static void mul(numtype p1[], const numtype p2[])
  {""" % (nvar,ndeg1,ndeg2)
    for k in range(len(terms_res)-1,-1,-1):
        pk = ["p1[%i]*p2[%i]" % (p[1],p[2]) for p in prods if p[0] == k]
        print "    p1[%i] =" % k," + ".join(pk),";"
    print """  }
};"""
    return True


def polymul_instance(nvar, ndeg1, ndeg2, maxterms = -1):
    terms1 = exponents(nvar,ndeg1)
    terms2 = exponents(nvar,ndeg2)
    terms_res = exponents(nvar,ndeg1+ndeg2)
    mon2 = monomial_exponents(nvar,ndeg2)
    prods = []
    mon_prods = []
    for (i1,t1) in enumerate(terms1):
        for (i2,t2) in enumerate(terms2):
            prodt = expadd(t1,t2)
            prods.append((terms_res.index(prodt),i1,i2))
    if maxterms > -1 and len(prods) > maxterms:
        return False
    for (i1,t1) in enumerate(terms1):
        for (i2,t2) in enumerate(mon2):
            prodt = expadd(t1,t2)
            mon_prods.append((terms_res.index(prodt),i1,i2))
    prods.sort()
    mon_prods.sort()
    print """template<class numtype>
struct polynomial_multiplier<numtype, %i, %i, %i>
{
  static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {""" % (nvar,ndeg1,ndeg2)
    for k in range(len(terms_res)):
        pk = ["p1[%i]*p2[%i]" % (p[1],p[2]) for p in prods if p[0] == k]
        print "    dst[%i] +=" % k," + ".join(pk),";"
    print """  }
  static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
  {"""
    for k in range(len(terms_res)):
        pk = ["p1[%i]*p2[%i]" % (p[1],p[2]) for p in prods if p[0] == k]
        print "    dst[%i] =" % k," + ".join(pk),";"
    print """  }"""
    print """   static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {""" 
    for k in range(len(terms_res)):
        pk = ["p1[%i]*m2[%i]" % (p[1],p[2]) for p in mon_prods if p[0] == k]
        if len(pk) > 0:
            print "    dst[%i] +=" % k," + ".join(pk),";"
    print """  }
  static void mul_monomial_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
  {"""
    for k in range(len(terms_res)):
        pk = ["p1[%i]*m2[%i]" % (p[1],p[2]) for p in mon_prods if p[0] == k]
        if len(pk) > 0:
            print "    dst[%i] =" % k," + ".join(pk),";"
    print """  }
  static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
  { 
  }
  static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
  {
  }    
};"""
    return True

print "// This file was automatically generated by gentab.py"
print 
for nvar in range(1,5):
    for ndeg in range(1,5):
        taylormul_instance(nvar,ndeg,ndeg,1000)
        print
        taylormul_inplace_instance(nvar,ndeg,ndeg,1000)
        print
        polymul_instance(nvar,ndeg,ndeg,1000)
        print
