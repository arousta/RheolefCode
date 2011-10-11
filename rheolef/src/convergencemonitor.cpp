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
#include "convergencemonitor.h"

//added by ali 11 Jun 2011
using rheolef::orheostream;

/**
 * @brief Default cpnstructor
 */
convergencemonitor::convergencemonitor()
:monitor(), current_line(0)
{
	
}

/**
 * @brief Constructor
 * @param headers a vector containing headers for the monitors
 */
convergencemonitor::convergencemonitor(vector<string> headers)
:current_line(0)
{
	vector<string>::iterator iter;
	for (iter = headers.begin(); iter < headers.end(); iter++)
	{
		monitor[*iter];
	}
}


/**
 * @brief Constructor
 * @param headers a vector containing headers for the monitors
 */
convergencemonitor::~convergencemonitor()
{
	
}

/**
 * @brief Copy constructor
 * @param rhs 
 */
convergencemonitor::convergencemonitor(const convergencemonitor& rhs)
:monitor(rhs.monitor),current_line(rhs.current_line)
{
	
}




/**
 * @brief Add a header to the monitor
 * @param header 
 */
void convergencemonitor::add_header(string header)
{
    monitor[header];
}






/**
 * @brief Add an entry line to monitor
 */
void convergencemonitor::add_line()
{
    current_line++;
	mapiterator iter;
	for( iter = monitor.begin(); iter !=  monitor.end(); iter++)
	{
		(iter->second).push_back(-1);
	}
}


/*!
    \fn convergencemonitor::add_entry(string entry, Float value)
	@brief Adds an entry to a monitor
 */
void convergencemonitor::add_entry(string entry, Float value)
{
    monitor[entry].push_back(value);
}


/*!
    \fn convergencemonitor::operator[](string)
	@brief Returns a reference to the corresponding monitor vector
 */
vector<Float> & convergencemonitor::operator[](string header)
{
    return monitor[header];
}


void write(const string & filename, const vector<string> & headers, const string & suffix = "conv")
{
	// Open iorheostream
    orheostream out_u(filename, suffix);

// 	for (vector<string>iterVstring = headers.begin(); iterVstring != headers.end();  ++iterVstring)
// 	{
// 		out_u << *iter_header << "\t" ;
// 		if ( monitors[*iter_header] != monitors.end() )
// 			size = max( size, *monitors[*iter_header].size() );
// 		else
// 			return;
// 	}
// 
// 	for ( int i = 0; i < size; i++ )
// 	{
// 		out_u << endl;
// 		for (iter_header = headers.begin(); iter_header != headers.end();  iter_header++)
// 		{
// 			cout	<< monitors[*iter_header]->at(i) << "\t";
// 		}
// 	}
	
    out_u.close();
	return;
}

