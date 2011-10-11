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
#ifndef SOLVER_PARTICLE_H
#define SOLVER_PARTICLE_H


//#include "util/XMLParser.h"
//#include "util/xmlfile.h"
#include "tensorfield.h"
#include "vectorfield.h"
#include "multifield.h"
#include "analytic_straightchannel.h"
#include "convergencemonitor.h"
#include "fluidmodel.h"

#include "rheolef.h"
#include "rheolef/uzawa_abtb.h"

//changed by ali 11 Jun 2011
#include "rheolef/uzawa_abtbc.h"

#include "rheolef/urm_abtbc.h"
#include "rheolef/eye.h"
#include "rheolef/trace.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <cmath>


#define EPS 1e-25

using namespace std;

//added by ali 11 Jun 2011
using rheolef::geo;
using rheolef::ssk;


/**
	@brief Solver for the wavychannel Problem.
	@author Andreas Putz
	This class implements the solver for the particulate flow problem.
	
	All the physical and computational parameters are stored in an xml file of the following form
	\verbatim
    <?xml version="1.0" encoding="UTF-8"?>
	<data>
		<parameters>
			<problem>ellipsoid</problem>
			<subproblem>gravity</subproblem>
			<workdir>.</workdir>
			<basename>ellipsoid</basename>

			<gridadaption>
				<gridadaptionsteps>1</gridadaptionsteps>
				<interactive> true | false </interactive>
				<criteria>gamma|gammadot|viscosity</criteria>
				<max_vertices>number</max_vertices>
				<bamg_options>
					<hmin>0.1</hmin>
					<hratio>1</hratio>
					<hcoef>1</hcoef>
				</bamg_options>
			</gridadaption>
	
			<FluidProperties>
				<Re>0</Re>
				<Bn>1</Bn>
				<bodyforce>0.</bodyforce>
			</FluidProperties>
	
			<ParticleProperties>
				<coefficient>-1</coefficient>
				<velocityguess>1</velocityguess>
				<FE_LagrangeMultSpace>P2</FE_LagrangeMultSpace>
			</ParticleProperties>
	
			<SolverViscoplastic>
				<type>Uzawa</type>
				<Model>Bingham</Model>
				<iterations>500</iterations>
				<augment>1e+7</augment>
				<tolerance>1e-5</tolerance>
				<FE_VelocitySpace>P2</FE_VelocitySpace>
				<FE_PressureSpace>P1</FE_PressureSpace>
				<FE_StressSpace>P1d</FE_StressSpace>
			</SolverViscoplastic>
			

			<StokesSolver>
				<type>default (Uzawa) | urm_abtbc </type>
				<iterations>50></iterations>
				<tolerance>1e-16</tolerance>
				<augment>1e+7</augment>
			</StokesSolver>
	
			<mesh>
				<geom>symxy</geom>
				<generator>bamg</generator>
				<area>1</area>
				<ratio>0.5</ratio>
				<psurface>0.1</psurface>
				<ymax>10</ymax>
				<xmax>10</xmax>
				<hmax>0.01</hmax>
				<subdivisions_sphere>100</subdivisions_top>
			</mesh>
		</parameters>
	</data>
	\endverbatim
	
	Parts of the data file:
 * - Viscoplastic Solver: 
 * 	\verbatim
		<SolverViscoplastic>
			<type>Uzawa</type>
			<iterations>1000</iterations>
			<augment>5e1</augment>
			<tolerance>1e-7</tolerance>
			<FE_VelocitySpace>P2</FE_VelocitySpace>
			<FE_PressureSpace>P1</FE_PressureSpace>
			<FE_StressSpace>P1d</FE_StressSpace>
		</SolverViscoplastic>
 * 	\endverbatim
 * 
 * 
 * 
	At the the parent node has to be set by the calling function "//parameters"
	@todo
		- Allocate space for the forms dynamically
		- Make the call more flexible
			- either set to "//parameters[last()]" to always select the last node
			- or allow to choose
		- Implement restart
		- Write each adaption step into the parameter.xml file with 
			\verbatim
			<parameters iteradapt=n>
				...
			</parameters>
			\endverbatim
		  or use this format for all the last childs ?
		- Write a decent convergence monitor class
			- use a dynamic matrix
			- vector for headers
			- file output
			- commandline output
			- convergence criteria
*/
class solver_particle{
private:
//	XMLParser*  xmlparam;
	//parameters params;
	int _timestep;
	Float _timereal;
	int _adaptstep;
	string _basename; ///< @brief Contains the filename with adaptive + time label
	string _geoname;  ///< @brief Contains the filename of the geometry without extension. 
	string _dataname; ///< @brief Contains the filename of the data without extension.

	geo omega_h;
	/** @name FE Spaces */
    //@{
	space Vh; ///< @brief Velocity space (vector space).
	space Vh0;///< @brief Velocity space (0 component).
	space Vh1;///< @brief Velocity space (1 component).
	space Qh; ///< @brief Pressure space.
	space Th; ///< @brief Stress, strain rate space (tensor space).
	space LPh;///< @brief Lagrange multiplier space for the particle
	space Bh; ///< @brief Combined Lagrange multiplier space ( pressure + particle )
	//@}
	/** @name FE Fields and Forms */
    //@{
	/**
	* @brief Map container for the finite element fields used repeatedly in the program.
 	*
	* Abandon to store every field.
	* -	u			velocity field (vector) P2
	* - p			pressure (scalar)		P1
	* - lambda 		ALE multiplier 			P2
	* - f			force field (rhs of the Stokes problem)
	* - gammadot	rate of shear tensor calculated from u
	* - gamma		relaxed rate of shear tensor
	* - T			Lagrange multiplier associated with the stress field.
	*/
	map<string, multifield> FEfields;
	typedef map<string, multifield>::iterator FEfieldsIterator;

	/**
	* @brief Map container for the finite element forms used repeatedly in the program.
 	*
	* Abandon to store every matrix.
	* Martices stored depend on the algorithm:
	* - Uzawa:
	*	- ar : Modified Laplace operator for the Stokes solver
			- Stokes Solver: default 
				- Weak Laplace operator a_h = 1.0*r_Bi*form(Vh,Vh,"2D_D"):
				  \f$ a(u,v):= 2 \int_\Omega D(u):D(v) dV \f$
				- ar = a_h + r_St*trans(b_h)*b_h
	*	- bh : grad, div operator -1*form(Vh,Qh,"div")
	*	- 2Du : form(Vh,Th,"2D")
	*	- DIV : -1*trans(2Du)
	*/
	map<string, form> FEforms;
	
	/**
	* @brief Stores the preconditioner for the Stokes problem
	*
	*/
	ssk<Float> fact;
	//@}

	fluidmodel _fluidmodel; ///< @brief Stores the fluid model

	/**
	* @brief Stores monitor vectors. Different lengths possible:
	*
	* - ellipsoid problem:
	*	"ParticleVelocity","a(u,u)","j(u)",["L(u)"]
	*/
	convergencemonitor monitors;

	/** 
	* @brief Contains the model parameters
	*
	* Possible values are:
	* - For Bingham const. law:
	*	- Bn Bingham number
	*/
	map<string, Float> _modelparams;

	/** 
	* @brief Contains counters to be appended to the _basename
	*/
	map<string, Float> _filenamecounters;
public:
//    solver_particle(XMLParser & xml_params);

    //new ctor by ali 11 Jun 2011
    solver_particle( std::string const& filename );

    ~solver_particle();
	void solve();
	///@name FE Methods
	//@{
	int FE_initialize();
	void FE_initialize_geo();
    int FE_initialize_spaces();
    int FE_initialize_fields();
    int FE_initialize_forms();
	//@}
	///@name Viscoplastic Solver: Uzawa Algorithms
	//@{
	void ALG_Uzawa_test(Float Lambda = 0.01);
    int ALG_Uzawa();
	//@}

	///@name Viscoplastic Solver: Regularisation Methods
	//@{
	Float VPREG_Lambda(int r, Float LambdaIn,Float LambdaTarget,string method = "times0.1");
	void VPREG_homotopy();
	//@}

	///@name Stokes Solver
	//@{
	int STOKES_Subproblem();
	int STOKES_Uzawa(field &);
	int Stokes_Piccard(bool updateprecond = true, Float = 0.01);
	int STOKES_urm_abtbc(field &, form &);
	field STOKES_rhs(field&,float);
	//@}
	void adapt_grid();
	///@name Set, Get, IO methods
	//@{
	void set_basename();
	void write_all();
	void write_fields();
	void write_fields(const string &);
	void write_geo();
	void write_geo(const string &);
	void write_stats(const string &);
    void plot_fields(string);
	void show_formstats(const form & testform, string ident="", string name="", bool plots = false);
	void show_spacestats(const space &, string, string);
    void analytic_solution(Float amplitude = 0, Float Bn = 0, Float ywall = 1,Float dPdL = 1);
	form mk_modtrace(string);
    
	//@}
};

#endif
