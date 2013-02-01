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

#ifndef STDMK_MESH_H
#define STDMK_MESH_H

// Standard C++ includes
#include <string>
#include <cmath>

// Rheolef Includes
#include "rheolef.h" 
#include "rheolef/uzawa_abtb.h"

#include "SimpleXML.h"

//added by ali 11 Jun 2011
using rheolef::Float;



using namespace std;
const Float pi=3.141592653589793238512808959406186204433;
/**
	@author Andreas Putz
	@brief The generator of the initial meshes.

	This class constructs the initial meshes using an external Grid
	Generator (BAMG or GRUMMP)

	All the physical and computational parameters are stored in an xml file of the following form
 * for the wavychannel problem
	\verbatim
    <?xml version="1.0" encoding="UTF-8"?>
	<data>
		<parameters>
			<problem>wavychannel</problem>
			<workdir>.</workdir>
			<basename>test</basename>
				
			<mesh>
				<geom>symy</geom>
				<generator>bamg</generator>
				<length>20</length>
				<height>1</height>
				<amplitude>-0.25</amplitude>
				<subdivisions_top>100</subdivisions_top>
				<bamg_options>
					<hmin>0.1</hmin>
					<hratio>1</hratio>
					<hcoef>1</hcoef>
				</bamg_options>
			</mesh>
		</parameters>
	</data>
	\endverbatim
	
	For the ellipsoid problem the following datafile is needed:
	\verbatim
	<?xml version="1.0" encoding="UTF-8"?>
	<data>
		<parameters>
			<problem>ellipsoid</problem>
			<subproblem>gravity</subproblem>
			<workdir>.</workdir>
			<basename>ellipsoid</basename>
	
			<mesh>
				<geom>symxy</geom>
				<generator>bamg</generator>
				<area>1</area>
				<ratio>0.5</ratio>
				<arclength>1.</arclength>
				<ymax>10</ymax>
				<xmax>10</xmax>
				<hmax>0.01</hmax>
				<subdivisions_sphere>100</subdivisions_top>
			</mesh>
		</parameters>
	</data>
	\endverbatim 

	@todo
	- Remove the grid visualisation from the generator methods and move 
	  them to a extra method.	
*/
class mk_mesh{
private:
	//gridparameters parameters;
    TiXmlDocument* xmldoc;
    TiXmlElement const* params;
    template< class T, int depth >
    void get_from_xml( path_t (&path)[depth], T*const x );

public:
	// Constructors Destructors
        mk_mesh( char const* fname );
        ~mk_mesh();

	//added by ali 11 Jun 2011
        int bamg_makemesh( );
        int bamg_wavychannel( );
        int bamg_ellipsoid( );
        //ali 12 Nov
        void bamg_squarechannel( );
        void bamg_trichannel();
        //ali 25 Jan 2013
        void bamg_bubblec();
};



template< class T, int depth >
inline void
mk_mesh::get_from_xml( path_t (&path)[depth], T *const x )
{get_value( path, params, x );}


#endif
