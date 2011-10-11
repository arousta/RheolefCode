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

#include "mk_mesh.h"


/**
   @brief Creates the standard meshes
   @param argc Number of input parameters
   @param argv Input Parameters
   @return error code
   
   \todo 
		- Handle multiple entries in input files
		- Call a general grid_construction method (not wavychannel)
 */
int zmain(int argc, char**argv)
{

  mk_mesh mesh( argv[1] );
  mesh.bamg_makemesh( );
  cout << "done" <<  endl;
	

  // Define Variables:
  //gridparameters parameters_grid;
  //mk_mesh mesh;

  //Create Testmesh
  //mesh.bamg_wavychannel(parameters_grid.wavychannel);
	
  return 1;
}
