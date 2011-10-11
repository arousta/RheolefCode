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
#ifndef VECTORFIELD_H
#define VECTORFIELD_H


#include "rheolef.h"
#include "rheolef/uzawa_abtb.h"

#include <string>
#include <iostream>

//added by ali 11 jun 2011
using rheolef::field;
using rheolef::space;
using rheolef::Float;

using namespace std;


/**
	@author Andreas Putz
	@brief Added support for tensors (Added support for stress tensors).
	This class adds a bit of support for tensors in the normal field class.
	@todo - Hand over tha mass form as a reference ? 
		  - Assignment operator needed
*/
class vectorfield : public field{
public:
    vectorfield();
	vectorfield(const space &V);
    ~vectorfield();
	vectorfield& operator=(const vectorfield&);
	vectorfield& operator=(const field&);
    Float norm(string, Float);
	Float dist(const vectorfield&,string, Float);
	void plot(string="NULL");
	///@name Analytic Solutions
	//@{
	void analytic_straightchannel(Float=0, Float=1, Float=1);
	void analytic_wavy_asymptotic(Float, Float, Float);
	//@}
    //field dotdot(const tensorfield &, const tensorfield &);

};

#endif
