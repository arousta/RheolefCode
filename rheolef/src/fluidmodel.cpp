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
#include "fluidmodel.h"

//added by ali 11 Jun 2011
using rheolef::sqr;


// fluidmodel::fluidmodel()
// {
// 	params["Bn"] = 0;
// 	_model = "Bingham";
// }


/**
 * @brief Default constructor
 */
fluidmodel::fluidmodel()
:_model("Simple"),_paramfield(NULL)
{
	params["Bn"] = 0;
}


/**
 * @brief Main constructor for the fluidmodel. Assigns a pointer to the reference field
 * @param inputfield Pointer to a reference field
 * @param model Simple, Bercovier, Papanastasiou
 */
fluidmodel::fluidmodel(multifield* inputfield, string model, const map<string, Float> & parameters)
:params(parameters),_model(model),_paramfield(inputfield)
{
}




fluidmodel::~fluidmodel()
{
}

/**
 * @brief Copy constructor
 * @param rhs Class to copy
 */
fluidmodel::fluidmodel(const fluidmodel& rhs)
{
	_model = rhs._model;
	_paramfield = rhs._paramfield;
	params = rhs.params;
}

/**
 * @brief Assignment operator
 * @param rhs Reference to another fluidmodel
 * @return 
 */
fluidmodel& fluidmodel::operator=(const fluidmodel& rhs)
{
	if (this != &rhs)
	{
		_model = rhs._model;
		_paramfield = rhs._paramfield;
		params = rhs.params;
	}
	return *this;
}


/*!
    \fn fluidmodel::set_params(const map<string, Float> parameters )
 */
void fluidmodel::set_params(const map<string, Float> & parameters )
{
   params = parameters;
}

void fluidmodel::set_modeltype(string model )
{
	_model = model;
}


/**
* @brief Returns a reference to the currently stored parameters
* @return reference to parameter map.
 */
 map< string, Float >& fluidmodel::getref_parameters()
{
	return params;
}


float fluidmodel::operator() (const point & x)
{
	Float eval;
	Float result;
	_paramfield->evaluate(x,eval);
	if(_model == "Simple")
	{	
		result = (1+ params["Bn"]/(params["Delta"]+eval) );
	}
	else if (_model == "Bercovier")
	{
		result = (1+ params["Bn"]/sqrt(sqr(params["Delta"])+sqr(eval)));
	}
	else if (_model == "Papanastasiou")
	{
		if (eval < 1e-25)
			result = 1+params["Bn"]/params["Delta"];
		else
			result = 1 + params["Bn"]*(1-exp(-eval/params["Delta"]))/eval;
	}
	else
	{
		cout << "Incorrect model" << endl;
	}
	return result;
}

/**

 * @brief Calculates the viscosity field corresponding to the model
 * @param[in,out] visc viscosity field
 * @param[in] rate of shear tensor field 
 */
void fluidmodel::viscosity(multifield & visc, const multifield & paramfield)
{
	//_paramfield = &paramfield;
// 	visc = interpolate(visc.get_space(),*this);
// 	for (multifield::size_type i=0; i<visc.n_component(); i++)
// 	{
// 		if (i==0)
// 			visc[0] = interpolate(paramfield->get_space(),*this);
// 		else
// 			visc[i] = visc[0];
// 	}
}
