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

#include "analytic_straightchannel.h"

//added by ali 11 Jun 2011
using rheolef::sqr;
using std::abs;
#include <cassert>

/**
 * @brief Default Constructor
 * @param Bn Bingham number
 * @param dpdx Pressure gradient \f$ \frac{\partial p}{\partial x_0} \f$
 * @param ywall coordinate of the wall
 */
analytic_straight_velocity::analytic_straight_velocity(Float Bn, Float dpdx, Float ywall) :
		_Bn(Bn), _dpdx(dpdx), _ywall(ywall)
{
}


/**
 * @brief This defines the functionpart of the whole thing. Needed for rheolef::interpolate
 *
 *	Let \f$ x = (x_0,x_1)\f$	

	\f[
		u(x) =  (u_0(x_1), 0)
	\f]
	\f{eqnarray}
		\begin{cases}
			\frac{f}{2}
			\left[
				|y_w|-\frac{Bn}{|f|}
			\right]^2
			&
			x_1 \in [\pm Bn / |f|] \\
			\frac{f}{2}
			\left[ 
				\left( |y_w| - \frac{Bn}{|f|} \right )^2-
				\left( |y|   - \frac{Bn}{|f|} \right )^2
			\right]
			&
			|x_1| >  Bn / |f|
		\end{cases}
	\f}
 * @param x point position
 * @return analytic velocity
 */
point analytic_straight_velocity ::operator() (const point & x)
{
	point u;
	Float yc = _Bn / abs(_dpdx);
	if (abs(x[1]) <= yc)
	{
		u[0]=(_dpdx/2)*sqr(abs(_ywall)-yc);
	}
	else
	{
		u[0]=(_dpdx/2)*(sqr(abs(_ywall)-yc)-sqr(abs(x[1])-yc));
	}
	return u;
}




//added by ali 24 Jul 2011
Float analytic_straight_velocity2::operator() (const point & x)
{
        point u;
        assert(x[1]<=1.);
        if (abs(x[1]) <= _ys)
        {
                u[0]=_Bn/(2.*_ys)*sqr(1.-_ys);
        }
        else
        {
                u[0]=_Bn/(2.*_ys)*(sqr(1.-_ys)-sqr(abs(x[1])-_ys));
        }
        return u[0];
}





/**
 * @brief Default Constructor
 * @param Bn Bingham number
 * @param dpdx Pressure gradient \f$ \frac{\partial p}{\partial x_0} \f$
 * @param ywall coordinate of the wall
 */
analytic_straight_shearrate::analytic_straight_shearrate(Float Bn, Float dpdx, Float ywall) :
		_Bn(Bn), _dpdx(dpdx), _ywall(ywall)
{
}

/**
 * @brief Definition of the abstract function operator.
 *
 * @param x Spatial position
 * @return Tensor containing the shearrate.
 */
tensor analytic_straight_shearrate ::operator() (const point & x)
{
	tensor gamma;
	Float yc = _Bn / abs(_dpdx);
	if (abs(_dpdx) > yc)
	{
		if (abs(x[1]) <= yc)
		{
			gamma=0;
		}
		else
		{
			gamma(0,0)=0;
			gamma(0,1)=-_dpdx*(x[1]-yc);
			gamma(1,0)=gamma(0,1);
			gamma(1,1)=0;
		}
	}
		return gamma;
}


/**
 * @brief Default Constructor
 * @param Bn Bingham number
 * @param dpdx Pressure gradient \f$ \frac{\partial p}{\partial x_0} \f$
 * @param ywall coordinate of the wall
 */
analytic_straight_stress::analytic_straight_stress(Float Bn, Float dpdx, Float ywall) :
		_Bn(Bn), _dpdx(dpdx), _ywall(ywall)
{
}

/**
 * @brief Definition of the abstract function operator.
 *
 * @param x Spatial position
 * @return Tensor containing the shearrate.
 */
tensor analytic_straight_stress ::operator() (const point & x)
{
	tensor T;
	
	T(0,0)=0;
	T(0,1)=-_dpdx*x[1];
	T(1,0)=T(0,1);
	T(1,1)=0;

	return T;
}
