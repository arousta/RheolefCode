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

#include "solver_particle.h"

// Command line parser

using namespace std;

/**
   @brief Solves the viscoplastic problem for ellipsoid particles
   @param argc Number of input parameters
   @param argv Input Parameters
   @return error code
   
   \todo 
		- Handle multiple entries in input files
 */
int main_d(int argc, char**argv)
{

// Wrap everything in a try block.  Do this every time, 
// because exceptions will be thrown for problems. 
	// Define the command line object.
//	CmdLine cmd("Solve the ellipsoid particle problem", ' ', "0.1");

	// Define the input file arguments
//	UnlabeledValueArg<string>  inputxml( "xmlfile", "Input XML file",true, "none","Input file (xml)"  );
//        cmd.add( inputxml );

	// Parse the commandline
//	cmd.parse( argc, argv );

	// Initialise the XML Class
//	string filename = inputxml.getValue();

//	XMLParser xml_parameters(filename);
//		cout << "Call Particle Solver class .." << endl;
//		solver_particle solver(filename); // previously: solver(xml_parameters);
//		solver.solve();
//		cout << "done" <<  endl;
	
  return 1;
}
