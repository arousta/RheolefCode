/***************************************************************************
 *   Copyright (C) 2007 by Andreas Putz   *
 *   putza@math.ubc.ca   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "solver_wavychannel.h"

#include <iostream>
#include <cstdlib>

using namespace std;

void read_xml( char const* fname );

/**
   @brief Solves the viscoplastic problem for wavy channels
   @param argc Number of input parameters
   @param argv Input Parameters
   @return error code
   
   \todo 
		- Handle multiple entries in input files
 */
int main(int argc, char**argv)
{

//  print_xml_doc( argc, argv );

  solver_wavychannel solver( argc, argv );
  solver.solve();

  return 1;
}
