/*
 * adaptive_strategy.h
 *
 *  Created on: Jul 7, 2011
 *      Author: ali
 */

#ifndef ADAPTIVE_STRATEGY_H_
#define ADAPTIVE_STRATEGY_H_

#include <vector>
#include <cstdlib>
#include <cassert>


struct adaptiveStrategy {

  double gamma_tol( size_t i )
  {
    check_range(i);
    return tol[i];
  }

  double hcoef( size_t i )
  {
    check_range(i);
    return coef[i];
  }

  size_t max_iterations( size_t i )
  {
    check_range(i);
    return maxit[i];
  }

  adaptiveStrategy()
  {
    double low = 1e-3;
    double med = 1e-4;
    double hig = 1e-6;
    double c[] = {      10,      9,      8,      7,      7,      7,      7};
    double t[] = {      low,    low,    med,    med,    hig,    hig,    hig};
    size_t m[] = {      1000,   1000,   1500,   2000,   2500,   3000,   4000};
    int nstep = 7;
    coef.assign (c, c+nstep);
    tol.assign( t, t+nstep);
    maxit.assign(m, m+nstep);
  }

private:
  std::vector<double> coef;
  std::vector<double> tol;
  std::vector<size_t> maxit;

  void check_range( size_t i )
  {assert( i<coef.size() );}
};



#endif /* ADAPTIVE_STRATEGY_H_ */
