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
#include "vectorfield.h"

//added by ali 12 Jun 2011
using namespace rheolef;

vectorfield::vectorfield()
{
}

/**
 * @brief Constructor taking the space as parameter. Zero preset
 * @param V space
 */
vectorfield::vectorfield(const space &V)
:field(V,0.)
{
}

vectorfield::~vectorfield()
{
}


/**
 * @brief Assignment operator for vectorfields
 * @param that reference to vectorfield to copy
 * @return reference to this
 */
vectorfield& vectorfield::operator=(const vectorfield& that)
{
	if (&that != this)
	{
		field::operator=(that);
	}
	return *this;
}


/**
 * @brief Assign a field to a vectorfield
 * @param that reference to vectorfield to copy
 * @return reference to this
 * @TODO check if the field is in fact vectorial
 */
vectorfield& vectorfield::operator=(const field& that)
{
	if (&that != this)
	{
		field::operator=(that);
	}
	return *this;
}


/**
 * \fn vectorfield::dot(const tensor &, const tensor &)
 * @brief Calculate the scalar inner procuct
 * @param cectorfield a
 * @param vectorfield b
 * @return a.b
 */
// field tensorfield::dot(const vectorfield &a__h, const vectorfield &b__h)
// {
// 	
// }


/**
 * @brief Calculate norms of the vectorfield
 * Calculates different norms of the vectorfield
 * Type has to be:
 * 		- Linfinity : Calculates \f$ \| u \|_\infty \f$ (mass_h is dummy)
 *		- L2		: \f$ \| u \|_{L2} \f$
 * @param form mass form
 * @param type type of the norm
 * @return Float norm of a
 */
Float vectorfield::norm(string type, Float scale = 1)
{
    if (type == "Linfinity"){
		return (*this).max_abs();
	}
	else if (type=="L2"){
	  const space& space_h = (*this).get_space();
	  form mass_h(space_h,space_h,"mass");
	  return sqrt(dot((*this),mass_h*(*this)));
	}
	else if (type=="flowrate"){
	  const space& space_h = (*this).get_space();
	  form mass_h(space_h,space_h,"mass");
	  field tt__h(space_h);
	  tt__h[0]=1;
	  tt__h[1]=0;
	  return scale*(dot(tt__h,mass_h*(*this)));
	}
	return -1.;
}


/**
 * @brief Calculate the distance between the current vectorfield and another vectorfield that.
 * Type has to be:
 * 		- Linfinity : Calculates \f$ \| u \|_\infty \f$ (mass_h is dummy)
 *		- L2		: \f$ \| u \|_{L2} \f$
 *		- flowrate  :
 *
 * @param that vectorfield
 * @param type string Linfinity/L2/flowrate
 * @param scale Float default = 1
 * @return \f$ \|this-that\| \f$
 */
Float vectorfield::dist(const vectorfield& that, string type, Float scale = 1)
{
	if (type == "Linfinity"){
		return (*this-that).max_abs();
	}
	else if (type=="L2"){
	  space space_h = (*this).get_space();
	  form mass_h(space_h,space_h,"mass");
	  return sqrt(dot((*this-that),mass_h*(*this-that)));
	}
	else if (type=="flowrate"){
	  space space_h = (*this).get_space();
	  form mass_h(space_h,space_h,"mass");
	  field tt__h(space_h);
	  tt__h[0]=1;
	  tt__h[1]=0;
	  return scale*(dot(tt__h,mass_h*(*this-that)));
	}
	return -1.;
}

void vectorfield::plot(string title)
{
	bool whilebreak=true;	
	int choice;
	do
	{
		cout<< "\t|===========================================================|" << endl
			<< "\t| Graphical Output : " << title << endl
			<< "\t| MENU:  0 exit        100 options" << endl
			<< "\t|===========================================================|" << endl
			<< "\t|  1 velocity vector field                                  |" << endl;
		for (size_type i=0; i<(*this).dimension(); i++)
		{
		cout<< "\t|  " << i+2 << " Component " << i << endl; 
		}
		cout<< "\t|-----------------------------------------------------------|" << endl;
		cin >> choice;

		switch (choice) {
		case 1:  { cout << velocity << plotmtv << *this; break; }
 		case 2:  
		  { 
			field uh_1 = (*this)[0]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
 		case 3:  
		  { 
			field uh_1 = (*this)[1]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 4:  
		  { 
			field uh_1 = (*this)[2]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 100:{	
					do{
					  cerr	<< endl <<"\t|************************************************************" << endl;
					  cerr 	<<"\t| Options for the graphical output "
							<< endl;
					  cerr	<< endl <<"\t|************************************************************" << endl;
					  cerr  << "\t| 1 Show FE Grid  |  2 Hide FE Grid  |" << endl;
					  cin >> choice ;
					  switch (choice) {
					    case  1 :{cout << grid; break;}
					    case  2 :{cout << nogrid; break;}
					    default :{choice = 0;}
					  }
				 	}
					while (choice);
					break;
				}
		default :{whilebreak = false;}
		}
	}while(whilebreak);
}



/**
 * @brief Calculates the analytic velocity solution for a straight channel
 * @param Float Bn Bingham number
 * @param Float deltaP pressure drop 
 * @param Float ywall coordinate of the wall
 */
void vectorfield::analytic_straightchannel(Float Bn, Float deltaP, Float ywall)
{
	//*this=interpolate(this->get_space(), analytic_straight_velocity(Bn,10,1));
}
