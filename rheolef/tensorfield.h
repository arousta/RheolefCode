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
#ifndef TENSORFIELD_H
#define TENSORFIELD_H


#include "rheolef.h"
#include "rheolef/uzawa_abtb.h"

#include <string>
#include <iostream>

#include "vectorfield.h"

//added by ali 11 Jun 2011
using rheolef::Form;

using namespace std;


/**
	@author Andreas Putz
	@brief Added support for tensors (Added support for stress tensors).
	This class adds a bit of support for tensors in the normal field class. 
	@TODO 
		- Default operators needed
			- Copy constructor
			- write clone fct
*/
class tensorfield : public field{
private:
	
public:
	///@name Constructors, Destructors, Assignment
	//@{
    tensorfield();
	tensorfield(const space &V);
	tensorfield(const tensorfield&);
    ~tensorfield();
	tensorfield& operator=(const tensorfield&);
	tensorfield& operator=(const field&);
	field& getfield();
	void u2shearrate(const vectorfield&, const form&, const form&);
	//@}
	///@name Constitutive Laws
	//@{
    tensorfield gammadot_to_stress(string constlaw, Float Bn);
	tensorfield stress_to_gammadot(string constlaw, Float Bn, Float scale);
	void set_stress_to_gammadot(string constlaw, Float Bn,tensorfield&, Float scale);
	//@}
	///@name Inner and outer products:
	//@{
	field dotdot(const tensorfield &);
    field secinv(const tensorfield &a__h);
	field secinv();
	//@}
	///@name IO
	//@{
	void plot(string="NULL");
	//@}
};

#endif
