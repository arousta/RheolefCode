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
#include "tensorfield.h"

//added by ali 12 Jun 2011
using namespace rheolef;


tensorfield::tensorfield()
:field()
{
}

/**
 * @brief Constructor taking the space as parameter. Zero preset
 * @param V space
 */
tensorfield::tensorfield(const space &V)
:field(V,0.)
{
}

// tensorfield::tensorfield(const tensorfield& that)
// {
// }

tensorfield::~tensorfield()
{
}



/**
 * @brief Assignment operator for tensorfields
 * @param that reference to tensorfield to copy
 * @return reference to this
 */
tensorfield& tensorfield::operator=(const tensorfield& that)
{
	if (&that != this)
	{
		field::operator=(that);
	}
	return *this;
}


/**
 * @brief Assign a field to a tensorfield
 * @param that reference to tensorfield to copy
 * @return reference to this
 * @TODO check if the field is in fact tensorial
 */
tensorfield& tensorfield::operator=(const field& that)
{
	if (&that != this)
	{
		field::operator=(that);
	}
	return *this;
}


/**
 * Calculates the rate of shear tensor from a given velocity
 * @param  const vectorfield& u__h	velocity field
 * @param  const form& inv_mass_h	inverse mass form
 * @param  const form& shearrate_h	shearrate form
 */
void tensorfield::u2shearrate(const vectorfield& u__h, const form& inv_mass_h, const form& shearrate_h)
{
	*this = (1.0*inv_mass_h*(shearrate_h*u__h));
}



/**
 * @brief Returns the reference to the parent class field.
 *
 * @return Reference to field
 */
field& tensorfield::getfield()
{
	field& that = field::operator=(*this);
	return that;
}

/**
 * \fn tensorfield::dotdot(const tensor &)
 * @brief Calculate the tensor inner procuct of a tensor b with the current tensor.
 * @param tensor b
 * @return (*this):b
 */
field tensorfield::dotdot(const tensorfield &b__h)
{
	if (((*this).dimension() == 2)&&(b__h.dimension() == 2))
	{
		return (*this)(0,0)*b__h(0,0)+2*(*this)(0,1)*b__h(0,1)+(*this)(1,1)*b__h(1,1);
	}
    
	else if (((*this).dimension() == 3)&&(b__h.dimension() == 3))
	{
		return (*this)(0,0)*b__h(0,0)+(*this)(1,1)*b__h(1,1)+(*this)(2,2)*b__h(2,2)
			   + 2*(*this)(0,1)*b__h(0,1)+(*this)(0,2)*b__h(0,2)
			   + 2*(*this)(1,2)*b__h(1,2);
	}
	else
	{
		//cerr << "Quantity ins not a tensorial field" << endl;
	}
	return 0.*(*this);
}

/**
    \fn tensorfield::secinv(const tensorfield &a__h)
 */
field tensorfield::secinv(const tensorfield &a__h)
{
 	return sqrt(0.5*dotdot(a__h));
}


/**
    \fn tensorfield::secinv()
	@brief Returns the second invariant of the tensor
 */
field tensorfield::secinv()
{
 	return secinv(*this);
}

/**
	@brief Computes the stress tensor from a given rate of shear tensor using a basic material law.
	@param tensorfield rate of shear tensor
	@param string	   Constitutive Law (Bingham)
	@param Float	   Bingham number
	Allowed values for string:
		-	Bingham
	@todo Hershley Bulkley
	@todo Newtonian
	@todo Return pointer and not the fiels -> faster
 */
tensorfield tensorfield::gammadot_to_stress(string constlaw, Float Bn)
{
    Float konst=0.;
	size_type i,j = 0;
	field::size_type complength = (*this).get_space().n_unknown_component(0);
	field second_invariant__h = secinv();
	tensorfield stress__h((*this).get_space());

	for (i=0; i<(*this).dimension(); i++)
	{
		if (constlaw == "Bingham")
		{
			for (j=i*complength;j<complength+i*complength;j++)
			{
				if (1. * second_invariant__h.u(j) <= Bn)
				{
					stress__h.u(j)=0;
				}
				else
				{
					konst= (1.+Bn/second_invariant__h.u(j));
					stress__h.u(j) = konst*(*this).u(j);
				}
			}
		}
	}
	return stress__h;
}



/**
 * @brief Computes the stress tensor from a given rate of shear tensor using a Bingham law.
 *
 * The Bingham law is slightly modified to allow the use for the augmented Lagrangian formulation.
 * @param constlaw 
 * @param Bn 
 * @param a 
 * @return rate of shear tensor
 * @todo 
 *	- Newtonian and HB law
 *  - Return pointer and not the fields -> faster
 */
tensorfield tensorfield::stress_to_gammadot(string constlaw, Float Bn, Float scale = 1)
{	
	//field secinv = secinv(stress__h);
	Float konst=0.;
	size_type i,j = 0;
	field::size_type complength = (*this).get_space().n_unknown_component(0);
	field second_invariant__h = secinv();
	tensorfield strainrate__h((*this).get_space());

	for (i=0; i<(*this).dimension(); i++)
	{
		if (constlaw == "Bingham")
		{
			for (j=i*complength;j<complength+i*complength;j++)
			{
				if (1. * second_invariant__h.u(j) <= Bn)
				{
					strainrate__h.u(j)=0;
				}
				else
				{
					konst= (1.-Bn/second_invariant__h.u(j));
					strainrate__h.u(j) = (konst/scale)*(*this).u(j);
				}
			}
		}
	}
	return strainrate__h;
}

/**
 * @brief Sets the tensorfield according to a scaled material law

 * - Bingham law
	 \f{align*}
		\gamma = 
		\begin{cases}
			0 
			& \text{if $\|T\|_{II}\leq Bn$} \\
			\left[
				1-\frac{Bn}{\|T\|_{II}}
			\right]
			\frac{T}{s}
		\end{cases}
	\f}
 * @param constlaw default: Bingham
 * @param Bn Bingham number
 * @param stress stress tensorfield \f$ T \f$
 * @param scale scale	\f$ s \f$
 */
void tensorfield::set_stress_to_gammadot(string constlaw, Float Bn,tensorfield& stress, Float scale)
{
	Float konst=0.;
	size_type i,j = 0;
	field::size_type complength = (stress).get_space().n_unknown_component(0);
	field second_invariant__h = stress.secinv();
	vector<field::size_type> lengthofcomps;
	size_type n_components = (stress).n_component();
	

	for (i=0; i<n_components; i++)
	{	
		lengthofcomps.push_back((stress).get_space().n_unknown_component(i));
	}

	if (constlaw == "Bingham")
	{
		for (j=0;j<lengthofcomps[0];j++)
		{
			if (1. * second_invariant__h.u(j) <= Bn)
			{
				int k=0;
				for(i=0; i<n_components;i++)
				{
					(*this).u(j+k)=0;
					k +=lengthofcomps[i];
				}
			}
			else
			{
				konst= (1.-Bn/second_invariant__h.u(j))/scale;
				int k = 0;
				for(i=0; i<n_components;i++)
				{
					(*this).u(j+k) = konst*stress.u(j+k);	
					k +=lengthofcomps[i];
				}
			}
		}
	}
	else
	{
		cerr << "No constitutive law given" << endl;
	}
}


void tensorfield::plot(string title)
{
	//cout << velocity << plotmtv << *this;
	bool whilebreak=true;	
	int choice;
	do
	{
		cout<< "\t|===========================================================|" << endl
			<< "\t| Graphical Output : " << title << endl
			<< "\t| MENU:  0 exit        100 options" << endl
			<< "\t|===========================================================|" << endl;
		for (size_type i=0; i<(*this).n_component(); i++)
		{
		cout<< "\t|  " << i+1 << " Component " << i << endl; 
		}
		cout<< "\t|  " << 9 << " Second invariant" << endl;
		cout<< "\t|-----------------------------------------------------------|" << endl;
		cin >> choice;

		switch (choice) {
 		case 1:  
		  { 
			field uh_1 = (*this)[0]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
 		case 2:  
		  { 
			field uh_1 = (*this)[1]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 3:  
		  { 
			field uh_1 = (*this)[2]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 4:  
		  { 
			field uh_1 = (*this)[3]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 5:  
		  { 
			field uh_1 = (*this)[4]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 6:  
		  { 
			field uh_1 = (*this)[5]; 
			cout << rheo << clean << fill;
			cout << plotmtv << uh_1;
			break;
		  }
		case 9:
		  {
			field uh_1 = (*this).secinv();
			cout << rheo <<fill << plotmtv << uh_1;
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
