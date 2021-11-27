#include <iostream>
#include "espic_math.h"

using namespace ESPIC;

namespace ESPIC{

Real newton_raphson(const std::function< Real(Real) >& f, 
                    const std::function< Real(Real) >& df)
{
  Real a = 1.0;
  Real rsd = f(a)/df(a);
  int it;
  bool converged = fabs(rsd) < rtol_nr;

  for (it = 0; !converged && it < maxit_nr; ++it) {
    a -= rsd;
    rsd = f(a)/df(a);
    converged = fabs(rsd) < rtol_nr;
  }

  if (converged) {
    std::cout << "Newton-Raphson: iter/maxit = " << it << "/" << maxit_nr
      << ": a = " << a << ".\n";
  }
  else {
    std::cout << "Newton-Raphson: failed to converge after " 
      << maxit_nr << " iterations.\n";
  }

  return a;
}

// define and initialize global seed for random number generator
Bigint Random::seed = 0;


}
