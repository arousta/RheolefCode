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
#ifndef SOLVER_FLOWONSET_H
#define SOLVER_FLOWONSET_H


#include "tensorfield.h"
#include "vectorfield.h"
#include "multifield.h"
#include "analytic_straightchannel.h"
#include "convergencemonitor.h"
#include "fluidmodel.h"

#include "rheolef.h"
#include "rheolef/uzawa_abtb.h"

//added by ali 12 Jun 2011
#include "rheolef/uzawa_abtbc.h"

#include "rheolef/urm_abtbc.h"
#include "rheolef/eye.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <cmath>

//added by ali 14 Jun 2011
#include "SimpleXML.h"
class TiXmlElement;

//added by ali 21 Oct 2011
#include "Restart.h"
#include "FieldVisualization.h"

#define EPS 1e-25

using namespace std;




/**
	@brief Solver for the wavychannel Problem.
	@author Andreas Putz

	This class implements the solver for the wavychannel problem.
	
	All the physical and computational parameters are stored in an xml file of the following form
	\verbatim
    <?xml version="1.0" encoding="UTF-8"?>
	<data>
		<parameters>
			<problem>wavychannel</problem>
			<subproblem>flowelementsrate | pressuredrop </subproblem>
			<workdir>.</workdir>
			<basename>test</basename>

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
			
			<flowrate>
				<dsteps>10</dsteps>
				<maxsteps>50</maxsteps>
				<target>1</target>
				<dPdL>2</dPdL>
				<tolerance>1e-5</tolerance>
			</flowrate>
			
			<pressuredrop>
				<dPdL>2</dPdL>
			</pressuredrop>
	
			<FluidProperties>
				<Re>0</Re>
				<Bn>1</Bn>
				<density>1000</density>
			</FluidProperties>
	
			<ParticleProperties>
				<density>2000</density>
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
			
			<SolverViscoplastic>
				<type>Uzawa</type>
				<Model>Simple</Model>
				<regparam_start>1e-2</regparam_start>
				<regparam_end>1e-19</regparam_end>
				<regparam_descend>times0.1</regparam_descend>
				<maxregsteps>10</maxregsteps>
				<iterations>500</iterations>
				<augment>1e+7</augment>
				<tolerance>1e-5</tolerance>
				<FE_VelocitySpace>P2</FE_VelocitySpace>
				<FE_PressureSpace>P1</FE_PressureSpace>
				<FE_StressSpace>P1d</FE_StressSpace>
			</SolverViscoplastic>			


			<SolverViscoplastic>
				<type>Regularised</type>
				<subtype>Piccard</subtype>
				<iterations>500</iterations>
				<Model>Simple | Bercovier | Papanastasiou </Model>
				<regparam_start>1e-2</regparam_start>
				<regparam_end>1e-19</regparam_end>
				<regparam_descend>times0.1</regparam_descend>
				<maxregsteps>10</maxregsteps>
				<nonlinmaxiter>100</nonlinmaxiter>
				<nonlinminiter>10</nonlinminiter>
				<nonlintol>1e-7</nonlintol>
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
				<geom>symy</geom>
				<generator>bamg</generator>
				<length>20</length>
				<height>1</height>
				<amplitude>-0.25</amplitude>
				<hmax>0.1</hmax>
				<subdivisions_top>100</subdivisions_top>
			</mesh>
		</parameters>
	</data>
	\endverbatim
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
//added by ali 12 Jun 2011
using rheolef::Float;
using rheolef::geo;
using rheolef::space;
using rheolef::ssk;

class solver_flowonset{
private:
        //added by ali 14 Jun 2011
        TiXmlElement const* parameters;
        TiXmlDocument* xmldoc;
        char const* xml_filename;
        void load_xml_file( );
        template< class T, int depth >
        void get_from_xml( path_t (&path)[depth], T*const x );
        template< int depth >
        bool path_exist( path_t (&path)[depth] );
        template< int depth >
        string get_string( path_t (&path)[depth] );
        //ali Oct 2011
        const char* restart_file;
        void restart();
        int uzawa_iter_restart;
        bool is_a_restart;
        Restart _restart;
        void save_monitors( std::ofstream& );
        void load_monitors( std::ifstream& );
        void write_restart_files( const int niter );
        void load_fields_for_postprocessing( const string& fieldfile );
        void viz1( );
        Postprocessing_FieldsOverLine horizn_output( const string& base, const string& title, const Float y );
        void put_fields_on_line( SegmentedLine& sline, Postprocessing_FieldsOverLine& r );
        //wn: width of narrower part of channel, ww: width of plug in wider part of channel
        bool is_plug_broken( Float& wn, Float& ww );
        // p & T on symmetry behave like ~xf(y), compute f(y) with a least square method
        void plot_fp_fT();
        //ali Feb 2012
        void dump_matlab( const char* fname );

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
	space Qh; ///< @brief Pressure space.
	space Th; ///< @brief Stress, strain rate space (tensor space).
	//@}
	/** @name FE Fields and Forms */
    //@{
	/**
	* @brief Map container for the finite element fields used repeatedly in the program.
 	*
	* Abandon to store every field.
	* -	u	velocity field (vector) P2
	* - p	pressure (scalar)		P1
	* - f	force field (rhs of the Stokes problem)
	* - gammadot	rate of shear tensor calculated from u
	* - gamma	relaxed rate of shear tensor
	* - T		Lagrange multiplier associated with the stress field.
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
				- bh : grad, div operator -1*form(Vh,Qh,"div")
				- ar = a_h + r_St*trans(b_h)*b_h
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
	* - flowrate problem:
	*	-  "dPdL" sequence of pressure drops
	*	-  "flowrate" corresponding pressure drops
	* - pressuredrop problem:
	*	- "dPdL"
	*/
	convergencemonitor monitors;

	/** 
	* @brief Contains the model parameters
	*
	* Possible values are:
	* - For Bingham const. law:
	*	- Bn Bingham number
	* - For Reg. Models (Simple, Bercovier, Pap)
	*	- Bn Bingham number
	*	- Delta reg. param
	*/
	map<string, Float> _modelparams;

	/** 
	* @brief Contains counters to be appended to the _basename
	*/
	map<string, Float> _filenamecounters;
public:

	//added by ali 12 Jun 2011
//    solver_flowonset(XMLParser & xml_params);
//	solver_flowonset( char const* filename );
	//ali 26 Oct 2011
	solver_flowonset( int argc, char** argv );

	~solver_flowonset();
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
	int STOKES_flowrate();
	int STOKES_Uzawa(field &);
	int Stokes_Piccard(bool updateprecond = true, Float = 0.01);
	int STOKES_urm_abtbc(field &, form &);
	field STOKES_rhs(field&,Float,bool,Float);
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
    void plot_fields(string);
    void analytic_solution(Float amplitude = 0, Float Bn = 0, Float ywall = 1,Float dPdL = 1);
    
	//@}
};





template< class T, int depth >
inline void
solver_flowonset::get_from_xml( path_t (&path)[depth], T *const x )
{get_value( path, parameters, x );}

template< int depth >
inline bool
solver_flowonset::path_exist( path_t (&path)[depth] )
{return xml_path_exist( path, parameters );}


template< int depth >
inline string
solver_flowonset::get_string( path_t (&path)[depth] )
{
  string s;
  get_from_xml( path, &s );
  return s;
}

#endif
