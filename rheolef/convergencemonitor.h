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
#ifndef CONVERGENCEMONITOR_H
#define CONVERGENCEMONITOR_H


#include "rheolef.h"

#include <string>
#include <iostream>
#include <map>
#include <vector>

//added by ali 11 Jun 2011
using rheolef::Float;

using namespace std;


/**
	@author Andreas Putz
	@brief Adds monitors for convergence parameters
	
*/
class convergencemonitor{
private:
	typedef map<string,vector<Float> >::iterator mapiterator;
	map<string,vector<Float> > monitor;
	
public:
	///@name Constructors, Destructors, Assignment
	//@{
	convergencemonitor();
	convergencemonitor(vector<string> headers);
	convergencemonitor(const convergencemonitor& rhs);
	~convergencemonitor();
	//@}

	///@name Get, set, add methods
	//@{
	void add_header(string header);
    void add_line();
    void add_entry(string entry, Float value);
    vector<Float> & operator[](string);

	
	//@}


	///@name IO Methods
	//@{
	void write(const string & filename, const vector<string> & headers, const string & suffix = "conv");
	//@}
protected:
    int current_line;
};

#endif
