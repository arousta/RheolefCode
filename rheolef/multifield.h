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
#ifndef MULTIFIELD_H
#define MULTIFIELD_H


#include "rheolef.h"
#include "rheolef/uzawa_abtb.h"

#include <string>
#include <iostream>


//added by ali 11 Jun 2011
using rheolef::field;
using rheolef::Form;
using rheolef::Float;
using rheolef::space;



/**
	@author Andreas Putz
	@brief Added support for tensors (Added support for stress tensors).
	This class adds a bit of support for tensors in the normal field class. 
	@todo 
		- Default operators needed
			- Copy constructor
			- write clone fct
*/
class multifield : public field{
private:

public:
	///@name Constructors, Destructors, Assignment
	//@{
    multifield();
	multifield(const space &V,Float = 0);
	multifield(const multifield&);
	multifield(const field&);
    ~multifield();

	multifield& operator=(const multifield&);
	multifield& operator=(const field&);

	field& getfield();
	void u2shearrate(const multifield&, const form&, const form&);
	//@}
	///@name Constitutive Laws, Derivatives
	//@{
    friend void gammadot_to_stress(multifield&, const multifield&, std::string, Float);
	friend void secinvgammadot_to_viscosity(multifield&, const multifield&, std::string, Float, Float);
	friend void stress_to_gammadot(const multifield&, multifield&, std::string, Float, Float=1, Float=0.01);
	friend void rateofstrain(multifield&, const multifield&,const form&, const form&, Float);
	//@}
	///@name Inner and outer products:
	//@{

	friend Float norm(const field& that,std::string normtype, Float scale = 1);
	friend Float flowrate(const field& u__h, const field& f__h,const form& mass_h ,Float scale = 1);
	friend field dotdot(const field &, const field &);
    friend field secinv(const field& a__h, Float scale = 0.5 );

	//@}
	///@name IO
	//@{
	void plot(std::string header = "" ) const;
	//@}
private:
	void plot_scalar(std::string title = "") const;
	void plot_vector(std::string title = "") const;
	void plot_tensor(std::string title = "") const;
};

#endif
