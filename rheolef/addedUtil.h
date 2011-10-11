/*
 * addedUtil.h
 *
 *  Created on: Jun 18, 2011
 *      Author: ali
 */

#ifndef ADDEDUTIL_H_
#define ADDEDUTIL_H_

#include <iostream>
#include <cstdio>

#include "rheolef/form.h"
#include "rheolef/field.h"
#include "rheolef/space.h"

using rheolef::Form;
using rheolef::field;
using rheolef::space;


inline void
print_DOF_numbers( form const& f, char const* name )
{
  printf("form %s:\n",name);
  printf("\tnrow: %u\n", f.nrow());
  printf("\tncol: %u\n", f.ncol());

  printf("uu block:\n");
  printf("\tnrow: %u\n", f.uu.nrow());
  printf("\tncol: %u\n", f.uu.ncol());

  printf("bb block:\n");
  printf("\tnrow: %u\n", f.bb.nrow());
  printf("\tncol: %u\n", f.bb.ncol());
  fflush(stdout);
}

inline void
print_DOF_numbers( field const& f, char const* name )
{
  printf("field %s, nDOF: %u\n", name, f.n_unknown()+f.n_blocked());
  printf("\tsize: %u\n",f.size());
  fflush(stdout);
}

inline void
print_DOF_numbers( space const& s, char const* name )
{
  printf("space %s, size:%u \n", name,  s.size());
  fflush(stdout);
}

inline void
raise_and_exit( char const* function, char const* file, size_t line, char const* err )
{
  printf("%s:\n\t%s @ %u: error\n\t  %s\n", file, function, line, err);
  exit(-1);
}


#endif /* ADDEDUTIL_H_ */
