/*
 * addedUtil.cpp
 *
 *  Created on: May 14, 2012
 *      Author: ali
 */

#include <iostream>
#include <iomanip>

#include "MyUtil.h"
#include "addedUtil.h"
#include "rheolef/csr.h"


void print( std::ostream& s, const rheolef::form& m )
{
  assert(bool(s));
  s.flags(s.scientific);
  s.precision(3);
  for( size_t i=0; i<m.nrow(); ++i ){
      for( size_t j=0; j<m.ncol(); ++j ){
          rheolef::Float entry;
          Myutil::set_Nan(&entry);
          const size_t uuR = m.uu.nrow();
          const size_t uuC = m.uu.ncol();
          const bool in_uurow = i<m.uu.nrow();
          const bool in_uucol = i<m.uu.ncol();
          if( in_uurow ){
              if( in_uucol )
                entry = m.uu(i,j);
              else
                entry = m.ub(i,j-uuC);
          }
          else {
              const size_t I = i-uuR;
              if( in_uucol )
                entry = m.bu(I,j);
              else
                entry = m.bb(I,j-uuC);
          }
          s << std::setw(12);
          s << entry;
      }
      s << '\n';
  }
}
