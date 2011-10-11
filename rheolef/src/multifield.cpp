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
#include "multifield.h"

using namespace std;

//added by ali 11 Jun 2011
using namespace rheolef;


/*!
	@brief Default constructor (Calls field constructor)
 */
multifield::multifield()
:field()
{
}

/*!
	@brief Constructor (Calls field constructor)
 */
multifield::multifield(const space &V,Float value )
:field(V,value)
{
}


/*!
	@brief Copy Constructor (Calls field constructor)
*/
multifield::multifield(const multifield& rhs)
:field(rhs)
{
}

/*!
	@brief Copy Constructor (Calls field constructor)
*/
multifield::multifield(const field& rhs)
:field(rhs)
{
}

/**
 * @brief Assignment operator for multifields
 * @param rhs reference to tensorfield to copy
 * @return reference to this
 */
multifield& multifield::operator=(const multifield& rhs)
{
	if (&rhs != this)
	{
		field::operator=(rhs);
	}
	return *this;
}


/**
 * @brief Assign a field to a multifield
 * @param rhs reference to tensorfield to copy
 * @return reference to this
 */
multifield& multifield::operator=(const field& rhs)
{
	if (&rhs != this)
	{
		field::operator=(rhs);
	}
	return *this;
}



/**
 * @brief Returns the reference to the parent class field.
 *
 * @return Reference to field
 */
field& multifield::getfield()
{
	field& that = field::operator=(*this);
	return that;
}






/**
 * @brief Plot the selected field
 * @param header Title of the graph
 */
void multifield::plot(std::string header) const
{
	if (this->get_valued()=="scalar")
	{
		plot_scalar(header);
	}
	else if (this->get_valued()=="vector")
	{
		plot_vector(header);
	}
	else if (this->get_valued()=="tensor")
	{
		plot_tensor(header);
	}
}


/**
 * @brief Plot a scalar field
 * @param title 
 */
void multifield::plot_scalar(std::string title) const
{
	bool whilebreak=true;	
	int choice;
	do
	{
		cout<< "\t|===========================================================|" << endl
			<< "\t| Graphical Output : " << title << endl
			<< "\t| MENU:  0 exit        100 options" << endl
			<< "\t|===========================================================|" << endl
			<< "\t|  1 Plot scalar field                                      |" << endl
			<< "\t|-----------------------------------------------------------|" << endl;
		cin >> choice;

		switch (choice) 
		{
			case 1:
			{
				cout << rheo << clean << fill;
				cout << plotmtv << *this;
			}
			case 100:
			{	
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
 * @brief Plot a vectorial multifield.
 * @param title 
 */
void multifield::plot_vector(std::string title) const
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
 * @brief Plot tensorial field
 * @param title 
 */
void multifield::plot_tensor(std::string title) const
{
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
			field uh_1 = secinv(*this);
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


/**
 * @brief Destructor (empty at the moment)
*/
multifield::~multifield()
{
}

/** 
 * @brief Calculate the tensor inner procuct of two tensor lhs,rhs
 * @param lhs tensorfield on the left side
 * @param rhs tensorfield on the right side
 * @return A:B
	
*/
field dotdot(const field & lhs, const field & rhs)
{
	if ((lhs.get_valued()=="tensor") && (rhs.get_valued()=="tensor"))
	{
		if ((lhs.dimension() == 2)&&(rhs.dimension() == 2))
		{
			return lhs(0,0)*rhs(0,0)+2*lhs(0,1)*rhs(0,1)+lhs(1,1)*rhs(1,1);
		}
		
		else if ((lhs.dimension() == 3)&&(rhs.dimension() == 3))
		{
			return lhs(0,0)*rhs(0,0)+lhs(1,1)*rhs(1,1)+lhs(2,2)*rhs(2,2)
				+ 2*lhs(0,1)*rhs(0,1)+lhs(0,2)*rhs(0,2)
				+ 2*lhs(1,2)*rhs(1,2);
		}
		else
		{
			cerr << "Quantity ins not a tensorial field" << endl;
		}
	}
	return 0.*lhs;
}


/**
	@brief Calculates the norm of a given field or multifield

	@param that field/multifield of which the norm will be calculated
	@param normtype	declares which norm to be taken:
		- "Linf"	\f$ \|u \|_{L^\infty}^2:= max_{x\in \Omega} |u(x)| \f$
		- "L2"		\f$ \|u \|_{L^2}^2:=
						\begin{cases}
							\int_\Omega |u|^2 dV & \text{if the field is scalar} \\
							\int_\Omega |u \circ u| dV & \text{if the field is vectorial}\\
							\int_\Omega 0.5 T:T dV & \text{if the field is tensorial}
						\end{cases}
					\f$
	@param scale scaling factor for the norms (default = 1)
	@return Norm of the field

	@todo
		- Write norm for:
			- tensorial:  Linf
			- H1 norm for all
		- Add a class which saves the mass matrix, perhaps in the class.
*/
Float norm(const field& that,std::string normtype, Float scale)
{
	if (that.get_valued()=="scalar" && normtype =="Linf")
	{
		return that.max_abs();
	} 
	else if ( (that.get_valued()=="scalar" || that.get_valued()=="vector") && normtype=="L2")
	{
		const space& space_h = that.get_space();
		form mass_h(space_h,space_h,"mass");
		return sqrt(dot(that,mass_h*that));
	}
	else if ( (that.get_valued()=="tensor") && normtype=="L2")
	{
		field tmp = 0.5*dotdot(that,that);
		const space& space_h = tmp.get_space();
		form mass_h(space_h,space_h,"mass");
		return sqrt(dot(tmp,mass_h*tmp));
	}
	else
	{
		cout 	<< "***************************************************************" << endl
				<< "* Error: Norm not defined" << endl
				<< "***************************************************************" << endl;
	}
	return -1;
}

/**
	@brief Calculate the fowrate of a given velocity field

	@param[in] velocity field
	@param[in] direction field
	@param[in] scalar product matrix
	@param	scale
	@return \f$ \|T\|_{II} := \sqrt{s T_{ij}T_{ij}} \f$
	
*/
Float flowrate(const field& u__h, const field& f__h,const form& mass_h ,Float scale)
{
	return scale*(dot(f__h,mass_h*u__h));
}


/**
	@brief Calculate the second invariant \f$\|T\|_{II} := \sqrt{s T_{ii}T_{ij}} \f$ of a given tensorial field T, where scale s is set to 0.5 by default
	
	@param a__h tensorfield
	@param scale Scale s for the second invariant (default = 0.5)
	@return \f$ \|T\|_{II} := \sqrt{s T_{ij}T_{ij}} \f$
	
*/
field secinv(const field& a__h, Float scale)
{
	return sqrt(scale*dotdot(a__h,a__h));
}


/**
	@brief Computes the stress tensor from a given rate of shear tensor using a basic material law.
	
	@param[out] stress_h Stress tensor
	@param[in] gammadot__h 	rate of strain tensor
	@param constlaw	Constitutive Law (Bingham). 
	Allowed values for constlaw are:
		-	Bingham \f$ T(x) := 
						\begin{cases}
							0 & \text{if } \|\dot \gamma (u(x)) \|_{II} = 0 \\
							\left( 1+ \frac{Bn}{\|\dot \gamma (u(x)) \|_{II}} \right) \dot \gamma (u(x))
							& \text{if } \|\dot \gamma (u(x)) \|_{II} > 0 \\
						\end{cases}
					\f$
	@param Bn	Bingham number
	@todo 
	- Different constitutive laws:
		- Newtonian
		- Herschel Bulkley
	- Loop over all components
	- Make sure that the both fields have the same space.
	- Don't rely on compenents being equally long.
 */
void gammadot_to_stress(multifield& stress__h,const multifield& gammadot__h, std::string constlaw, Float Bn)
{
	Float konst=0.;
	field::size_type i,j = 0;
	field::size_type complength = gammadot__h.get_space().n_unknown_component(0);
	field second_invariant__h = secinv(gammadot__h);

	for (i=0; i<gammadot__h.dimension(); i++)
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
					stress__h.u(j) = konst*gammadot__h.u(j);
				}
			}
		}
	}
}



/**
 * @brief Calculates the viscosity of the regularised model w.r.t the second invariant of the shearrate.

 * @param viscosity__h reference to the viscosity field
 * @param sgammadot__h reference to the rate of shear field
 * @param constlaw Constitutive Law (Simple, Bercovier, Papanastasiou)
 * @param Bn Bingham number
 * @param reg regularisation parameter
 */
void secinvgammadot_to_viscosity(multifield& viscosity__h,const multifield& sgammadot__h, std::string constlaw, Float Bn, Float reg)
{
	multifield::size_type nunknowns = sgammadot__h.n_unknown(), i;
	cerr << reg << endl;
	if (nunknowns != viscosity__h.n_unknown() )
	{
		cout << " Different degrees of freedom within the two fields" << endl;
	}
	
	if (constlaw == "Simple")
	{
		for (i=0; i<nunknowns; i++)
		{
			viscosity__h.u(i) =  1+ Bn/( reg + sgammadot__h.u(i) );
		}	
	}
	

	
// 	Float konst=0.;
// 	field::size_type i,j = 0;
// 	field::size_type complength = sgammadot__h.get_space().n_unknown_component(0);
// 	
// 	for (i=0; i<sgammadot__h.dimension(); i++)
// 	{
// 		if (constlaw == "Simple")
// 		{
// 			for (j=i*complength;j<complength+i*complength;j++)
// 			{
// 				viscosity__h.u(j) = 1+ Bn/( reg + sgammadot__h.u(j) );
// 			}
// 		}
// 	}
}





/**
 * @brief Calculates the rate of strain tensor from a given stress tensor

 * @param[in] stress__h Stress tensor
 * @param[out] gammadot__h Rate of strain tensor
 * @param constlaw default: Bingham, possible values are
	- Bingham law
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
 * @param Bn Bingham number
 * @param stress stress tensorfield \f$ T \f$
 * @param scale Scale s. Important for the augmented Lagrange alg.
 * @param Lambda regularisation parameter
 *
 * @todo
	- Modify arithmetic
	- Make sure blocked parts are treated
	- Make sure that both spcaes are the same
 */
void stress_to_gammadot(const multifield& stress__h, multifield& gammadot__h, std::string constlaw, Float Bn, Float scale, Float Lambda )
{
	Float konst=0.;
	field::size_type i,j = 0;
	field::size_type complength = stress__h.get_space().n_unknown_component(0);
	field second_invariant__h = secinv(stress__h);
	vector<field::size_type> lengthofcomps;
	field::size_type n_components = (stress__h).n_component();

	for (i=0; i<n_components; i++)
	{	
		lengthofcomps.push_back((stress__h).get_space().n_unknown_component(i));
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
					gammadot__h.u(j+k)=0;
					k +=lengthofcomps[i];
				}
			}
			else
			{
				konst= (1.-Bn/second_invariant__h.u(j))/scale;
				int k = 0;
				for(i=0; i<n_components;i++)
				{
					gammadot__h.u(j+k) = konst*stress__h.u(j+k);	
					k +=lengthofcomps[i];
				}
			}
		}
	}
	else if (constlaw == "Simple")
	{
		for (j=0;j<lengthofcomps[0];j++)
		{
			Float konst0 = second_invariant__h.u(j) - Bn -scale*Lambda;
			
			konst= (konst0 + sqrt( 4*second_invariant__h.u(j)*Lambda*scale + sqr(konst0) ) )/
					(2*second_invariant__h.u(j)*scale);
			int k = 0;
			for(i=0; i<n_components;i++)
			{
				gammadot__h.u(j+k) = konst*stress__h.u(j+k);	
				k +=lengthofcomps[i];
			}
		}
	}
	else
	{
		cerr << "No constitutive law given" << endl;
	}
}




/**
 * @brief Calculates the rate of shear tensor from a given vectrofield

 * @param[out]  gammadot__h reference to the rate of shear field. Avoids multiple copying.
 * @param  u__h	velocity field
 * @param  inv_mass_h	inverse mass form
 * @param  shearrate_h	shearrate form
 * @param  scale scale to switch between $D(u)$ and rate of shear definition
 * @return rate of strain tensor
 */
void rateofstrain(multifield& gammadot__h, const multifield& u__h,const form& inv_mass_h, const form& shearrate_h, Float scale = 1)
{
	gammadot__h = scale*inv_mass_h*(shearrate_h*u__h);
}
