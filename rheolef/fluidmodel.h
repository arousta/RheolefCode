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
#ifndef FLUIDMODEL_H
#define FLUIDMODEL_H

#include "rheolef.h"
#include "vector"
#include "multifield.h"

#include <string>
#include <iostream>
#include <map>
#include <vector>

//added
using rheolef::Float;
using rheolef::point;


using namespace std;
/**
@brief Implements the fluid models

Implements the following fluid models:
- Bingham 
- Regularised Bingham models

	@author Andreas Putz <putza@math.ubc.ca>
*/
struct fluidmodel : unary_function< point, Float >
{
private:
	typedef map< string, Float > maptype;
	typedef map< string, Float > ::iterator iterparams;
	map<string, Float > params;
    string _model;
	multifield* _paramfield;
public:
	///@name Constructors, Destructors, Assignment
	//@{
    fluidmodel();
	fluidmodel(multifield *, string, const map<string, Float> &);
	fluidmodel(const fluidmodel&);
    ~fluidmodel();
	fluidmodel& operator=(const fluidmodel& );
	//@}
	float operator() (const point & x);

    void set_params(const map<string, Float> & parameters );
    void set_modeltype(string model="Simple" );
    maptype& getref_parameters();
    void viscosity(multifield & visc, const multifield & paramfield);

};

#endif
