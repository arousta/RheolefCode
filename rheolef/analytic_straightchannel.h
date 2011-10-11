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
#ifndef ANALYTIC_STRAIGHTCHANNEL_H
#define ANALYTIC_STRAIGHTCHANNEL_H

#include "vector"
#include "rheolef.h"

#include "tensorfield.h"
#include "vectorfield.h"

//added by ali 11 Jun 2011
using rheolef::point;
using rheolef::tensor;
#include <cmath>

/**
	* @author Andreas Putz <putza@math.ubc.ca>
	* @brief Calculates the analytic velocity to the pressure driven plain channel.
	* This can be used together with the interpolate function or rheolef.
	*
	* The expected geometry is a full or half channel. If a full channel is provided
	* it should range in [-ywall,ywall]
*/
struct analytic_straight_velocity : unary_function< point, point >
{
  public:
	analytic_straight_velocity(Float Bn = 0, Float dpdx = 1, Float ywall = 1);
	point operator() (const point & x);
  protected:
	Float _Bn;		///< @brief Bingham number
	Float _dpdx;	///< @brief pressure gradient \f$ \frac{\partial p}{\partial x_0} \f$.
	Float _ywall;	///< @brief Coordinate of the wall.
};


//added by ali 24 Jul 2011
//this is a modified version
struct analytic_straight_velocity2 : unary_function< point, Float >
{
  public:
        analytic_straight_velocity2( Float ys ):
          _Bn( 2./((std::pow(ys,3)+2.)/(3.*ys)-1.) ), _ys(ys)
        {
          cerr << "Computed Bn is: " << _Bn << '\n\n';
        }

        Float operator() (const point & x);
  protected:
        Float _Bn;              ///< @brief Bingham number
        Float _ys;
};

//added by ali 29 Jul 2011
struct newtonian_prof : unary_function< point, Float >
{
  Float operator()( const point& x ){
    return 1.5*( 1.-std::pow(x[1],2) );
  }
};


/**
	* @author Andreas Putz <putza@math.ubc.ca>
	* @brief Calculates the analytic shearrate to the pressure driven plain channel.
	* This can be used together with the interpolate function or rheolef.
	*
	* The expected geometry is a full or half channel. If a full channel is provided
	* it should range in [-ywall,ywall]
	* 
	* Let \f$ x = (x_0,x_1)\f$	be the position vector, and let

	\f[
		\dot \gamma(u) = \nabla u + (\nabla u)^t =  
		\begin{pmatrix} 
			\gamma_{0,0} & \gamma_{0,1} \\
			\gamma_{1,0} & \gamma_{1,1}
		\end{pmatrix}
	\f]
	be the symmertic rate of shear tensor. Then the analyic solution for the shearrate
	of a plane channel flow is
	\f{align}
	  \dot \gamma(u) = 
		\begin{cases}
			0
			& \text{if $\frac{\partial p}{\partial x_0} \leq  \frac{Bn}{|y_w|}$} \\
			\begin{cases}
				0
				& \text{if $|x_1| <  \frac{Bn}{|y_w|} $} \\
				\begin{pmatrix}
					0 & -\frac{\partial p}{\partial x_0} \left(x_1-\frac{Bn}{|y_w|} \right) \\
					-\frac{\partial p}{\partial x_0} \left( x_1-\frac{Bn}{|y_w|} \right) & 0
				\end{pmatrix}
				& \text{if $|x_1| \geq  \frac{Bn}{|y_w|} $}
			\end{cases}
			& \text{if $\frac{\partial p}{\partial x_0} >  \frac{Bn}{|y_w|}$}
		\end{cases}
	\f}
*/
struct analytic_straight_shearrate : unary_function< point, tensor >
{
  public:
	analytic_straight_shearrate(Float Bn = 0, Float dpdx = 1, Float ywall = 1);
	tensor operator() (const point & x);
  protected:
	Float _Bn;		///< @brief Bingham number
	Float _dpdx;	///< @brief pressure gradient \f$ \frac{\partial p}{\partial x_0} \f$.
	Float _ywall;	///< @brief Coordinate of the wall.
};

/**
	* @author Andreas Putz <putza@math.ubc.ca>
	* @brief Calculates the analytic stress tensor to the pressure driven plain channel.
	* This can be used together with the interpolate function or rheolef.
	*
	* The expected geometry is a full or half channel. If a full channel is provided
	* it should range in [-ywall,ywall]
	* 
	* Let \f$ x = (x_0,x_1)\f$	be the position vector, and let

	\f[
		\dot T = 
		\begin{pmatrix} 
			T_{0,0} & T_{0,1} \\
			T_{1,0} & T_{1,1}
		\end{pmatrix}
	\f]
	be the symmertic rate of stress tensor. Then the analyic solution for the shearrate
	of a plane channel flow is
	\f{align}
	  T = 
		\begin{pmatrix}
			0 & -\frac{\partial p}{\partial x_0} x_1 \\
			-\frac{\partial p}{\partial x_0} x_1 & 0
		\end{pmatrix}
	\f}
*/
struct analytic_straight_stress : unary_function< point, tensor >
{
  public:
	analytic_straight_stress(Float Bn = 0, Float dpdx = 1, Float ywall = 1);
	tensor operator() (const point & x);
  protected:
	Float _Bn;		///< @brief Bingham number
	Float _dpdx;	///< @brief pressure gradient \f$ \frac{\partial p}{\partial x_0} \f$.
	Float _ywall;	///< @brief Coordinate of the wall.
};
#endif
