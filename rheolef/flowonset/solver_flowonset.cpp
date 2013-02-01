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

/**
 * \file solver_flowonset.cpp
 * \author Andreas Putz
 * \date 01-01-2007
 * \brief Class implementation of the wavychannel problem for viscoplastic (Bingham) fluids
 *
 * This program implements the main solver class for the wavychannel problem.
*/


#include "solver_flowonset.h"
#include "SimpleXML.h"
#include "addedUtil.h"
//ali 7 Jul 2011
#include "adaptive_strategy.h"
#include "analytic_straightchannel.h"
//adaptiveStrategy adapt_strategy;
//int i_adapt = 0;
#include "FieldVisualization.h"
//ali 21 Oct 2011
#include "Restart.h"
#include "FieldVisualization.h"
#include <cstdlib>
#include <iomanip>
//ali 16 Jan 2012
#include "MyUtil.h"
#include <utility>
#include <cassert>
//ali 17 Feb 2012
#include <iostream>
#include <fstream>

using namespace rheolef;


void solver_flowonset::load_xml_file( )
{
  if( !(xmldoc->LoadFile(xml_filename)) ){
    cout << "Can not load the file " << xml_filename << endl;
    exit(1);
  }
  parameters = xmldoc->RootElement()->FirstChildElement("parameters");
}


/**
 * @brief Constructor initialising the xmlparam
 * @param xml_parameters
 */
solver_flowonset::solver_flowonset( int argc, char** argv )
:monitors()
,parameters(0) /*added by ali 16 Jun 2011*/
,xml_filename(argv[1])
,xmldoc( new TiXmlDocument )
,restart_file("restart_info") //ali 24 Oct 2011
,uzawa_iter_restart(0)
,is_a_restart(false)
{
    cout << "\t==========================================\n"
    	 << "\t|Constructor flowonset                   |\n"
	 << "\t==========================================" << endl;

    //added by ali 14 Jun 2011
    load_xml_file();
    set_basename();

    _timestep=0;
    _timereal=0.;
    _adaptstep=0;

    //ali 25 Oct 2011
    {
      path_t path[] = { "SolverViscoplastic", "restart_save_period" };
      double period;
      get_from_xml(path,&period);
      _restart.set_restart_saving_period(period);
    }
    if( Restart::restart_file_exists(restart_file) ){
        cout << "\t|\n"
             << "\t| Restart file found. Continueing from last save...\n"
             << "\t|@ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @" << endl;
        is_a_restart = true;
        restart();
        return;
    }

    FE_initialize();
	
	cout << "\t| - Setup pressure and flowrate ... \n";
	//added by ali 15 Jun 2011
	//xmlparam->get_childelement_string("subproblem")
	std::string str_subproblem;
	path_t subproblem_PT[] = {"subproblem"};
	get_from_xml( subproblem_PT, &str_subproblem );
	if (str_subproblem == "onset")
	{	
            cout << "\t|   Initialise onset problem\n";
            {
              Float P0;
              path_t flow_rate_dpdl_PT[] = {"onset","startP"};
              get_from_xml( flow_rate_dpdl_PT, &P0 );
              monitors["dPdL"].push_back(P0);
            }
	}

	cout 	<< "\t|Constructor flowonset completed         |\n"
		<< "\t==========================================" << endl << endl;
}


solver_flowonset::~solver_flowonset()
{
  //added by ali 16 Jun 2011
  delete xmldoc;
}


void solver_flowonset::FE_initialize_geo()
{
	if (_adaptstep == 0)
	{
            cout << "\t| - Initialize GEO ... \n";
            omega_h = geo(_geoname+".geo");
            cout << "\t|    *done\n"
                 << "\t-------------------------------------------" << endl;
	}
}


/**
 * \fn solver_flowonset::FE_initialize_spaces()
 * @return flag
 * @brief Initialize the FE spaces (order, block DOF)
 */
int solver_flowonset::FE_initialize_spaces()
{
    // Initialize Spaces
	//cout << plotmtv << omega_h;
    string velspace,pspace,tspace;
    cout << "\t| - Setup spaces ..." << endl;

    //added by ali 15 Jun 2011
    {
      path_t path[] = {"SolverViscoplastic","FE_VelocitySpace"};
      get_from_xml(path, &velspace);
    }
    {
      path_t path[] = {"SolverViscoplastic","FE_PressureSpace"};
      get_from_xml(path, &pspace);
    }
    {
      path_t path[] = {"SolverViscoplastic","FE_StressSpace"};
      get_from_xml(path, &tspace);
    }

    cout << "\t|   VelSpace = " << velspace  << ", PSpace = " << pspace << ", StressSpace = " << tspace << endl;

    Vh = space(omega_h,velspace,"vector");
    Qh = space(omega_h,pspace);
    Th = space(omega_h,tspace,"tensor");

    string geom;
    //added by ali 15 Jun 2011
    {
      path_t path[] = {"onset","geo"};
      get_from_xml(path, &geom);
//      geom = "full";
    }

    // Block for boundary conditions
    if (geom == "symy" || geom == "symxy")
    {
        Vh.block("top"); 	// Block top
        Vh[1].block("bottom"); 	// Symmetry     (uy = 0)
        Vh[1].block("left");	// Inflow       (uy = 0)
        Vh[1].block("right");	// Outflow      (uy = 0)
    }
    else if (geom == "full")
    {
        Vh.block("top"); 	// Block top
        Vh.block("bottom"); 	// Block bottom
        Vh[1].block("left");	// Inflow       (uy = 0)
        Vh[1].block("right");	// Outflow      (uy = 0)
    }
        //added by ali for testingggg
        //INFLOW
	// Block pressure space left and right:
	Qh.block("left");
	//if symxy is used don't block right added by ali for testing
	Qh.block("right");
	Qh.block("bubble");

        //INFLOW
	//For testinggg added by ali
//	Vh[0].block("left");

    cout << "\t|   Geometry type: " << geom << endl;
    cout << "\t|   *done\n" 
	 << "\t-------------------------------------------" << endl;
    return 1;
}


/*!
    \fn solver_flowonset::FE_initialize_fields()
 */
int solver_flowonset::FE_initialize_fields()
{
      cout << "\t| - Setup fields ..." << endl;
    string geom;
    string ViscoplSolver;
    //added by ali 15 Jun 2011
//#    xmlparam->get_childelement_stream("SolverViscoplastic/type") >>  ViscoplSolver;
    path_t solver_type_PT[] = {"SolverViscoplastic","type"};
    get_from_xml( solver_type_PT, &ViscoplSolver );

    if (ViscoplSolver == "Uzawa")
    {
          FEfields["u"] 	  = multifield(Vh,0);
          FEfields["f"]	          = multifield(Vh,0);
          FEfields["p"] 	  = multifield(Qh,0);
          FEfields["gamma"]       = multifield(Th,0);
          FEfields["gammadot"]    = multifield(Th,0);
          FEfields["T"] 	  = multifield(Th,0);
    }
    else if (ViscoplSolver == "Regularised")
    {
          FEfields["u"]              = multifield(Vh,0);
          FEfields["f"]	             = multifield(Vh,0);
          FEfields["p"] 	     = multifield(Qh,0);
          FEfields["gammadot"]       = multifield(Th,0);
          FEfields["secinv_gammadot"]= secinv(FEfields["gammadot"]);
          FEfields["viscosity"]      = multifield(FEfields["secinv_gammadot"].get_space(),1);
    }
    else
    {
          cout << "\t|    No valid method for the viscoplastic solver found" << endl;
    }

    //added by ali 24 Jul 2011 for INFLOW testing
    //INFLOW
    {
//    analytic_straight_velocity2 profile(0.223462071669235);
//    newtonian_prof profile;
//    domain const& inflow = omega_h["left"];
//    space Wh( omega_h, inflow, "P2" );
//    FEfields["u"][0]["left"] = interpolate( Wh, profile );
    }


    //added by ali 15 Jun 2011
//#	xmlparam->get_childelement_stream("mesh/geom") >> geom;
    path_t mesh_geo[] = {"mesh","geom"};
    get_from_xml(mesh_geo, &geom);

  // Impose unit force field
    //15 Jan 2013, added geom=="full", I guess we should do this for all geos?
    if (geom=="symy" || geom=="symxy" || geom=="full")
    {
	cout << "\t|    geom = " << geom << endl
	     << "\t|    Construct f = (1,0,0)^T" << endl;
        FEfields["f"][0] = 1;
	FEfields["mf"]   = FEforms["m"]*FEfields["f"];
    }

    return 1;
}


/*!
    \fn solver_flowonset::ALG_Uzawa()
	@brief Uzawa Algorithm for the augmented Lagrangian problem
	@todo Add description
 */
int solver_flowonset::ALG_Uzawa()
{
	bool FLAG_DEBUG_GRAPH = false;
    // ===========================================================
    // Other local parameters
    // ===========================================================
    int iter_uzawa=0;
    //added by ail 15 Jun 2011 types were float and I replaced with double
    Float dPdL_current = 1;	///< @todo Replace with flowrate from xml
    {
     path_t pressuredrop_dpdl_PT[] = {"pressuredrop","dPdL"};
     get_from_xml( pressuredrop_dpdl_PT, &dPdL_current );
    }


    Float dPdL_prev,flowrate_prev,flowrate_current,flowrate_target;
    multifield mf__h; 		// Temporary scalar product with force/pressure term

    // temporary field
    multifield tt__h;

    // -----------------------------------------------------------

    // ==========================================================
    // Set local parameters from the XML File
    // ==========================================================
    double r_Bi;
    Float Bn;
    int iter_uzawa_max,dflowrate;
    bool is_flowrate = false;
    string stokes_alg;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
//#    xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//#    xmlparam->get_childelement_stream("SolverViscoplastic/iterations") >> iter_uzawa_max;
//#    xmlparam->get_childelement_stream("StokesSolver/type") >> stokes_alg;
    //added by ali 15 Jun 2011
    {
      path_t path[] = {"SolverViscoplastic","augment"};
      get_from_xml(path, &r_Bi);
    }
    {
      path_t path[] = {"FluidProperties","Bn"};
      get_from_xml(path, &Bn);
    }
    {
      path_t path[] = {"SolverViscoplastic","iterations"};
      get_from_xml(path, &iter_uzawa_max);
    }
    {
      path_t path[] = {"StokesSolver","type"};
      get_from_xml(path, &stokes_alg);
    }

    //added by ali 15 Jun 2011
    //xmlparam->parse_childelement("subproblem")
    path_t subproblem_PT[] = {"subproblem"};
    if (path_exist(subproblem_PT))
    {
        string subproblem;
//#        xmlparam->get_data_stream() >> subproblem;
        get_from_xml( subproblem_PT, &subproblem );

        if (subproblem == "flowrate")
        {
            is_flowrate=true;
            //added by ali 15 Jun 2011
//#            xmlparam->get_childelement_stream("flowrate/dsteps") >> dflowrate;
//#            xmlparam->get_childelement_stream("flowrate/target") >> flowrate_target;
            {
              path_t path[] = {"flowrate","dsteps"};
              get_from_xml(path, &dflowrate);
            }
            {
              path_t path[] = {"flowrate","target"};
              get_from_xml(path, &flowrate_target);
            }
            cout << "\t|   Flowrate problem (steps: " << dflowrate
                 << ", target: "<< flowrate_target << ")" << endl;
        }
    }
    // -----------------------------------------------------------
//    analytic_solution(amplitude, Bn, ywall , dPdL_current);

    // ===========================================================
    // Main Uzawa Loop
    // ===========================================================

    for (iter_uzawa=0;iter_uzawa<iter_uzawa_max;iter_uzawa++)
    {
        cout << "\t|   * timestep = " << _timestep 
             << ", Adapt # = " << _adaptstep
             << ", Iteration # = " << iter_uzawa << endl;

        // =======================================================
        // Step 1 (Stokes Problem)
        // =======================================================
		cout << "\t|     - Step 1 (Stokes Problem)" << endl;

        field mrhs__h =	 STOKES_rhs(mf__h,r_Bi,is_flowrate,dPdL_current); // Assemble rhs
		if(false)
		{
			cout << "\t|       Plot mf__h" << endl;
			cout << mayavi << mf__h;
			cout << "\t|       Plot mrhs" << endl;
			cout << mayavi << mrhs__h;
		}
        // Use the Uzawa algorithm to solve the system. Use the inverse of A as preconditioner
        /// @todo Add selection of different solvers
        //if ( stokes_alg == "default" || stokes_alg == "Uzawa")
            STOKES_Uzawa(mrhs__h);
//            cout << velocity << gnuplot << FEfields["u"];
		//Calculate shearrate
        FEfields["gammadot"] = FEforms["u2shear"]*FEfields["u"];
		if(FLAG_DEBUG_GRAPH)
			{
				plot_fields("Step 1: Velocity and derived rate of shear");
			}


        // -------------------------------------------------------

        // =======================================================
        // Step 2 (Viscoplastic Problem)
        // =======================================================
		cout << "\t|     - Step 2 (Viscoplastic Problem)" << endl;
        // The use of the temprary stress is propably not optimal,
        // so it should be removed by implementing a conversion
        // method in the tensorfield class
        tt__h = FEfields["T"]+r_Bi*FEfields["gammadot"];
        //added by ali 21 Jun 2011 changed to FEfields["gamma"] first was FEfields["gammadot"]
		stress_to_gammadot( tt__h, FEfields["gamma"], "Bingham", Bn, 1+r_Bi);
		if(FLAG_DEBUG_GRAPH)
			{
				plot_fields("Step 2: new gammadot");
			}
        // -------------------------------------------------------

        // =======================================================
        // Step 3 (Update of Lagrange Mult)
        // =======================================================
		cout << "\t|     - Step 3 (Lagrange Mult)" << endl;
		FEfields["T"] = FEfields["T"] + 
			r_Bi*(FEfields["gammadot"]-FEfields["gamma"]);

        // -------------------------------------------------------

        // =======================================================
        // Step 4 (Flowrate step)
        // =======================================================
		
        if ( is_flowrate && (iter_uzawa % dflowrate == 0 ) && false)
        {
            cout << "\t|     - Step 4 (Flowrate step)" << endl;
            flowrate_current = norm(FEfields["u"],"flowrate",1);
            dPdL_current = dPdL_current
                             - (flowrate_current-flowrate_target)
                             *(dPdL_current-dPdL_prev)
                             /(flowrate_current-flowrate_prev);
        }

        // -------------------------------------------------------
    }
    return 1;
}


/*!
    \fn solver_flowonset::ALG_Uzawa_test()
	@brief Uzawa Algorithm for the augmented Lagrangian problem. Use the analytic solution.
	@todo Add description
 */
void solver_flowonset::ALG_Uzawa_test( Float Lambda )
{
	int FLAG_DEBUG_GRAPH = -1;
	cout <<"\t| - ALG_Uzawa" << endl;
    // ===========================================================
    // Other local parameters
    // ===========================================================
    int iter_uzawa=0;

    multifield tt__h; // temporary field
	

    // ==========================================================
    // Set local parameters from the XML File
    // ==========================================================
    Float r_Bi, Bn, tol, tolerance;
    int iter_uzawa_max;
    bool isflowrate = false;
    string stokes_alg;
    //added by ali 15 Jun 2011
    {
      path_t path[] = {"SolverViscoplastic","augment"};
      get_from_xml(path, &r_Bi);
    }
    {
      path_t path[] = {"FluidProperties","Bn"};
      get_from_xml(path, &Bn);
    }
    {
      path_t path[] = {"SolverViscoplastic","iterations"};
      get_from_xml(path, &iter_uzawa_max);
    }
    //by ali 7 Jul 2011
//    iter_uzawa_max = adapt_strategy.max_iterations(i_adapt);
    {
      path_t path[] = {"StokesSolver","type"};
      get_from_xml(path, &stokes_alg);
    }
    {
      path_t path[] = {"SolverViscoplastic","tolerance"};
      get_from_xml(path, &tolerance);
    }
    //ali 25 Oct 2011
    int nreport;
    {
      path_t path[] = {"SolverViscoplastic","nreport"};
      get_from_xml(path,&nreport);
    }
    //ali 19 Jan 2013
    Float Dp, onset_tol;
    {
      path_t flow_rate_dpdl_PT[] = {"onset","Pincrement"};
      get_from_xml( flow_rate_dpdl_PT, &Dp );
    }
    {
      path_t flow_rate_dpdl_PT[] = {"onset","tolerance"};
      get_from_xml( flow_rate_dpdl_PT, &onset_tol );
    }

    //added by ali 7 Jul 2011
//    tolerance = adapt_strategy.gamma_tol(i_adapt);

    //ali 24 Oct 2011
     if( is_a_restart )
     {
         iter_uzawa = uzawa_iter_restart;
         //set to false for later mesh adaption iterations
         is_a_restart = false;
     }
     tol = 1.e3;


    //Ali 17 Jan 2013
    bool resetLagMult(false);
    {
      std::string reset;
      path_t path[] = {"onset","resetLagMult"};
      get_from_xml(path, &reset);
      if(reset=="yes")
        resetLagMult = true;
    }
    bool enhanceLoop(false);
    bool flow_onset(false);

while (!flow_onset){

    if( resetLagMult ){
      FEfields["T"] =     0*FEfields["T"];
      FEfields["gamma"] = 0*FEfields["gamma"];
    }

    do
    {
        if(iter_uzawa!=0)
          cout << iter_uzawa << ",";
		multifield oldu(FEfields["u"]);
        // =======================================================
        // Step 1 (Stokes Problem)
        // =======================================================
		int StokesStatus = STOKES_flowrate();
		
		FEfields["gammadot"] = FEforms["u2shear"]*FEfields["u"];
		if(((iter_uzawa % FLAG_DEBUG_GRAPH) == 0 && FLAG_DEBUG_GRAPH > 0))
			{
				plot_fields("Step 1: Velocity and derived rate of shear");
			}
		
 
        // =======================================================
        // Step 2 (Viscoplastic Problem)
        // =======================================================
		tt__h = FEfields["T"] + r_Bi*FEfields["gammadot"];
		//added by ali 12/16 Jun 2011
		std::string str_SolverViscoplastic_Model;
		path_t solverVisc_model_PT[] = {"SolverViscoplastic","Model"};
		get_from_xml( solverVisc_model_PT, &str_SolverViscoplastic_Model );

		if (str_SolverViscoplastic_Model =="Bingham")
		{
//			cout << "Bingham Model" << endl;
			stress_to_gammadot( tt__h, FEfields["gamma"], "Bingham", Bn, 1.0+r_Bi );
		}
		else if (str_SolverViscoplastic_Model =="Simple")
		{
			cout << "Simple Model" << endl;
			stress_to_gammadot( tt__h, FEfields["gamma"],"Simple", Bn, 1.0+r_Bi, Lambda );
		}
		else
		{
			cerr << "No fluid model specified";
		}

		if((iter_uzawa % FLAG_DEBUG_GRAPH) == 0 && FLAG_DEBUG_GRAPH > 0)
			{
				plot_fields("Step 2: new gammadot");
			}
        

        // =======================================================
        // Step 3 (Update of Lagrange Mult)
        // =======================================================
		FEfields["T"] = FEfields["T"] + 
			r_Bi*(FEfields["gammadot"]-FEfields["gamma"]);

   
		if( (iter_uzawa % FLAG_DEBUG_GRAPH) == 0 && FLAG_DEBUG_GRAPH > 0)
			{
				plot_fields("Field Output after all Uzawa substeps");
			}
		// Convergence:
		//added by ali 6 Jul 2011
		if( iter_uzawa%nreport==0 )
		  {
		FEfields["secinv_gammadot"]=secinv(FEfields["gammadot"]);
		FEfields["secinv_gamma"]=secinv(FEfields["gamma"]);
		Float tolu = norm(multifield(FEfields["u"]-oldu),"L2");
		Float toldgamma = norm(multifield(FEfields["secinv_gamma"]-FEfields["secinv_gammadot"]),"L2");
		tol = max(tolu,toldgamma);

		cout  <<"\n\t\t    ||u-oldu||_L2           = " << tolu
		     << "\n\t\t    ||gamma-gammadot||_L2   = " << toldgamma << endl;
                if(tolerance<tol)
                  cout << "\tIterations #= ";
		  }
		iter_uzawa++;
    } while(iter_uzawa<=iter_uzawa_max && tol > tolerance);


    //ali 15 Jan 2013
    iter_uzawa = 0;
    FEfields["secinv_gamma"]=secinv(FEfields["gamma"]);
    const Float secnorm = rheolef::norm(FEfields["secinv_gamma"]);
    const bool zero_sec = (secnorm==0);
    Float curdPdL = monitors["dPdL"].back();
    if(enhanceLoop){
      Dp /= 2;
      if(Dp<=onset_tol)
        flow_onset = true;
    }
    if( zero_sec ){
        monitors["dPdL"].push_back(curdPdL+Dp);
    }
    else {
        if(!enhanceLoop){ // runs only once
          enhanceLoop = true;
          Dp /= 2;
        }
        monitors["dPdL"].push_back(curdPdL-Dp);
    }
    if(Dp==0)
      flow_onset = true;
    cout << "Setting dPdL: ========================( " << monitors["dPdL"].back() << " )" << endl;
}
}

/*!
    \fn solver_flowonset::FE_initialize_forms()
	@brief Initialise the FE operator forms
 */
int solver_flowonset::FE_initialize_forms()
{
    cout << "\t| - Setup forms ..." << endl;
    /////////////////////////////////////////////////////////////////////////
    // Setup the integral forms (new way)
    /////////////////////////////////////////////////////////////////////////
	
    string type;
    //added by ali 12 Jun 2011
//#	xmlparam->get_childelement_stream("SolverViscoplastic/type") >> type;
    //added by ali 15 Jun 2011
    {
      path_t path[] = {"SolverViscoplastic","type"};
      get_from_xml(path, &type);
    }


    if (type == "Uzawa")
    {
	cout << "\t|   Viscoplasitic Solver: " << type;

	form a_h(Vh, Vh, "2D_D");
	form shearrate_tensor_h(Vh,Th,"2D");		// $int_\Omega \bs \shearr(\mb u):T dv$
	form inv_mass(Th,Th,"inv_mass"); 		// Inverse mass operator for 2 tensors
	form mass_h(Th,Th,"mass");			// Tensorial mass form for L2 norm
	
	FEforms["b"] = form(Vh, Qh, "div");
	FEforms["b"] = -FEforms["b"]; 			// grad, divergence operator
	
	FEforms["m"] 		= form(Vh,Vh,"mass");	// Scalar product \int_Omega f \cdot v
	FEforms["DIV_ThVh"] = -1. * trans(shearrate_tensor_h);
	FEforms["u2shear"]	= inv_mass*shearrate_tensor_h;

	//added by ali 15 Jun 2011
	//xmlparam->parse_childelement("SolverViscoplastic/augment")
	path_t path[] = {"SolverViscoplastic","augment"};
    	if (path_exist(path))
        {
    	    double r_Bi;
//#         xmlparam->get_data_stream() >> r_Bi;
    	    get_from_xml(path, &r_Bi);

            a_h = 1.0*r_Bi*a_h;
            cout << ", r_Bi: " << r_Bi << endl;
        }

        //added by ali 15 Jun 2011
        //xmlparam->parse_childelement("StokesSolver")
    	path_t stokes_path[] = {"StokesSolver"};
        // Setup fields for the Stokes solver
        if (path_exist(stokes_path))
        {
            //added by ali 15 Jun 2011
            //xmlparam->parse_childelement("StokesSolver/type")
            path_t stokes_type_path[] = {"StokesSolver","type"};
            if (path_exist(stokes_type_path))
            {
//#		  xmlparam->get_data_stream() >> type;
                  get_from_xml(stokes_type_path, &type);
                  cout << "\t|   Stokes Solver: " << type;
                  if (type == "Uzawa" || type == "urm_abtbc")
                  {
                        //added by ali 15 Jun 2011
                        //xmlparam->parse_childelement("StokesSolver/augment")
                        path_t stokes_aug[] = {"StokesSolver","augment"};
                        if (path_exist(stokes_aug))
                        {
                            double r_St;
                            //added by ali 15 Jun 2011
//#		           xmlparam->get_data_stream() >> r_St;
                            get_from_xml(stokes_aug, &r_St);

                            FEforms["ar"] = a_h + r_St*trans(FEforms["b"])*FEforms["b"];
                            fact = ldlt(FEforms["ar"].uu);
                            cout << ", r_St: " << r_St << endl;
                        }
                  }
                  else
                  {
                       cerr << "No valid Stokes solver specifier" << endl;
                  }
            }//stokes solver type
        }//stokes solver
    }
    else if (type == "Regularised")
    {
        cout << "\t|   Viscoplasitic Solver: " << type;
        //FEforms["shear_T"] = form(Vh,Th,"2D");	// $int_\Omega \bs \shearr(\mb u):T dv$

        form shearrate_tensor_h(Vh,Th,"2D");	// $int_\Omega \bs \shearr(\mb u):T dv$
        form inv_mass(Th,Th,"inv_mass"); 	// Inverse mass operator for 2 tensors
// 	form mass_h(Th,Th,"mass");		// Tensorial mass form for L2 norm

        FEforms["b"] = form(Vh, Qh, "div");
        FEforms["b"] = -FEforms["b"]; 		// grad, divergence operator

        FEforms["m"] = form(Vh,Vh,"mass");	// Scalar product \int_Omega f \cdot v
        FEforms["u2shear"]	= inv_mass*shearrate_tensor_h;
    }
    else
    {
        cerr <<  "No vicoplastic solver specified" << endl;
         return -1;
    }

    cout << "\t|   *done\n"
         << "\t-------------------------------------------" << endl;
    return 1;
}



/*!
 * @brief Initialize the finite element spaces, forms and fields
 *
 */
int solver_flowonset::FE_initialize()
{
    cout << "\t******************************************\n"
         << "\t|FE initialize                           |\n"
         << "\t******************************************" << endl;
    FE_initialize_geo();
    FE_initialize_spaces();
    FE_initialize_forms();
    FE_initialize_fields();
    cout << "\t|FE initialized                          |\n"
         << "\t******************************************" << endl;
    return 1;
}




/**
 * @brief Solves the Stokes problem with the flowrate constraint
 *
 * @param  target target flowrate
 * @param  target_accuracy accuracy of the flowrate
 * @return 1 solver converged inside the step limit
 * @return -1 solver divergent
 * @todo implement me !!!
 */
int solver_flowonset::STOKES_flowrate()
{
	field tmpfield;
	//added by ali 15 Jun 2011
	std::string str_SolverViscoplastic_type;
	path_t solverVisc_type_PT[] = {"SolverViscoplastic","type"};
	get_from_xml( solverVisc_type_PT, &str_SolverViscoplastic_type );

	if (str_SolverViscoplastic_type=="Uzawa")
		{
			Float r_Bi;
			//added by ali 12 Jun 2011
			path_t solverVisc_aug_PT[] = {"SolverViscoplastic","augment"};
			get_from_xml( solverVisc_aug_PT, &r_Bi );

			tmpfield =  0.5*FEforms["DIV_ThVh"]*(FEfields["T"]-1*r_Bi*FEfields["gamma"]);
		}

		// Calculate the right hand side of the problem:
		field forcefield = *(monitors["dPdL"].end()-1)*FEfields["mf"];
		if (str_SolverViscoplastic_type == "Uzawa" )
			forcefield = forcefield + tmpfield;

		//added by ali 12/15 Jun 2011
		path_t stokes_type_PT[] = {"StokesSolver","type"};
		std::string str_StokesSolver_type;
		if (get_string(stokes_type_PT) == "Uzawa" )
			STOKES_Uzawa( forcefield );

	return 1;
}


/*!
 * @brief Solves the Stokes problem for a given right hand side.
 * (Stabilized version: cg_abtbc)
 *
 * \f[
	\begin{pmatrix}
		ar & b^T \\
		b  & C
	\end{pmatrix}
	\begin{pmatrix}
		u \\
		p 
	\end{pmatrix}
	=
	\begin{pmatrix}
		f \\
		g 
	\end{pmatrix}
   \f]
 * With a matrix preconditioner of the form
 * \f[
	\begin{pmatrix}
		ar^{-1} & \\
			    & I 
	\end{pmatrix}
   \f]
 * Solve a generalized Stokes problem of the form
 * @param rhs__h Finite element form of the rhs
 * @todo
	- Finish documentation
 */
int solver_flowonset::STOKES_urm_abtbc(field & rhs__h, form & c_h)
{
    double r_Stokes;
    int iter_Stokes_max;
    double tol_Stokes;
    double h;

    //added by ali 12/15 Jun 2011
//#    xmlparam->get_childelement_stream("StokesSolver/augment")    >> r_Stokes;
//#    xmlparam->get_childelement_stream("StokesSolver/iterations") >> iter_Stokes_max;
//#    xmlparam->get_childelement_stream("StokesSolver/tolerance")  >> tol_Stokes;
    {
      path_t path[] = {"StokesSolver","augment"};
      get_from_xml(path, &r_Stokes);
    }
    {
      path_t path[] = {"StokesSolver","iterations"};
      get_from_xml(path, &iter_Stokes_max);
    }
    {
      path_t path[] = {"StokesSolver","tolerance"};
      get_from_xml(path, &tol_Stokes);
    }

    h = omega_h.hmin();
	
    uzawa_abtbc (FEforms["ar"].uu, fact, FEforms["b"].uu, c_h.uu, 
				 FEfields["u"].u, FEfields["p"].u, rhs__h.u - (FEforms["ar"].ub*FEfields["u"].b), -(FEforms["b"].ub*FEfields["u"].b),
				 r_Stokes,iter_Stokes_max, tol_Stokes);
    return 1;
}

/*!
 * @brief Solves the Stokes problem for a given right hand side.
 *
 * Solve a generalized Stokes problem of the form
 * \f[
	\begin{pmatrix}
		ar & b^T \\
		b  & 0
	\end{pmatrix}
	\begin{pmatrix}
		u \\
		p 
	\end{pmatrix}
	=
	\begin{pmatrix}
		f \\
		g 
	\end{pmatrix}
   \f]
 * @param rhs__h Finite element form of the rhs
 * @todo
	- Finish documentation
 */
int solver_flowonset::STOKES_Uzawa(field & rhs__h)
{
    double r_Stokes;
    int iter_Stokes_max;
    double tol_Stokes;

    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("StokesSolver/augment")    >> r_Stokes;
//#    xmlparam->get_childelement_stream("StokesSolver/iterations") >> iter_Stokes_max;
//#    xmlparam->get_childelement_stream("StokesSolver/tolerance")  >> tol_Stokes;
    {
      path_t path[] = {"StokesSolver","augment"};
      get_from_xml(path, &r_Stokes);
    }
    {
      path_t path[] = {"StokesSolver","iterations"};
      get_from_xml(path, &iter_Stokes_max);
    }
    {
      path_t path[] = {"StokesSolver","tolerance"};
      get_from_xml(path, &tol_Stokes);
    }

// 	cout << "DOF u: unknowns:" << FEfields["u"].n_unknown()
// 		 << ", blocked: " << FEfields["u"].n_blocked() << endl;
// 	print_DOF_numbers(FEforms["ar"],"ar");

    uzawa_abtb (FEforms["ar"].uu, fact, FEforms["b"].uu,
				FEfields["u"].u, FEfields["p"].u, rhs__h.u - (FEforms["ar"].ub*FEfields["u"].b),  	-(FEforms["b"].ub*FEfields["u"].b), 
				r_Stokes, iter_Stokes_max, tol_Stokes);

    return 1;
}

/**
 * @brief Calculates the right hand side of the stoke problem
 *
 * @param mf__h
 * @param r_Bi
 * @param is_flowrate
 * @param dPdL
 * @return
 */
//added by ali 20 Jun 2011
//change the function prototype from
//field solver_flowonset::STOKES_rhs(field& mf__h, float r_Bi, bool is_flowrate = false, float dPdL=0)
field solver_flowonset::STOKES_rhs(field& mf__h, Float r_Bi, bool is_flowrate, Float dPdL)
{
	//cout << "\t|\t Construct rhs";
    field mdiv__h = FEforms["DIV_ThVh"]*(FEfields["T"]-1*r_Bi*FEfields["gamma"]);
//    cout << gnuplot << mdiv__h;
    if (is_flowrate==true)
    {
        FEfields["f"][0] = dPdL;
        mf__h = FEforms["m"]*FEfields["f"];
    }
	else
	{
		if (mf__h.dimension() <2)
		{
			FEfields["f"][0] = dPdL;
        	mf__h = FEforms["m"]*FEfields["f"];
		}
	}
    return 0.5*mdiv__h+mf__h;
}



/**
 * @brief Adapt the grid
 *
 * This method uses either the shear rate tensor of the Lagrange multiplier
 * to conduct (bamg: anisotropic) grid refinement.
 * Also reinterpolates the velocity
 * @todo Reread XML file
 */
void solver_flowonset::adapt_grid()
{
    string is_interactive;
    string criteria;
    Float Bn;
    field criteria__h;
    geo tmp_h;

    //added by ali 6 Jul 2011
    FEfields["secinv_gammadot"]=secinv(FEfields["gammadot"]);



	cout 	<< "\n**************************************************" << endl
        	<< "* Adaptive step                                  *" << endl;
	bool flag = false;
	do
	{
		flag = true;
		//added by ali 12/15 Jun 2011
//#		xmlparam->reparse_file();
//#		xmlparam->get_childelement_stream("gridadaption/interactive") >> is_interactive;
//#    	        xmlparam->get_childelement_stream("gridadaption/criteria") >> criteria;
		{
		  path_t path[] = {"gridadaption","interactive"};
		  get_from_xml(path, &is_interactive);
		}
		{
		  path_t path[] = {"gridadaption","criteria"};
		  get_from_xml(path, &criteria);
		}
		adapt_option_type options;

		{
		  path_t path[] = {"FluidProperties","Bn"};
		  get_from_xml(path, &Bn);
		}
		{
		  path_t path[] = {"gridadaption","max_vertices"};
		  get_from_xml(path, &options.n_vertices_max);
		}
		{
		  path_t path[] = {"gridadaption","bamg_options","hcoef"};
		  get_from_xml(path, &options.hcoef);
		}
		//added by ali 7 Jul 2011
//		options.hcoef = adapt_strategy.hcoef(i_adapt);
		{
		  path_t path[] = {"gridadaption","bamg_options","hmin"};
		  get_from_xml(path, &options.hmin);
		}
		{
		  path_t path[] = {"gridadaption","bamg_options","hmax"};
		  get_from_xml(path, &options.hmax);
		}
		{
		  path_t path[] = {"gridadaption","bamg_options","hratio"};
		  get_from_xml(path, &options.ratio);
		}


		// Select the adaption criteria
		if (criteria == "gamma")
		{
			field tmp = secinv(FEfields["gamma"]);
			criteria__h = sqrt(sqr(tmp)+Bn*tmp);
		}
		else if (criteria == "gammadot")
			criteria__h = sqrt(sqr(FEfields["secinv_gammadot"])+Bn*FEfields["secinv_gammadot"]);
		else if (criteria == "viscosity")
			criteria__h = FEfields["viscosity"];
		else
			cout << "No adaption criteria selected" << endl;
		// insert break condition

		// Create temporary geometry:
		tmp_h   = omega_h;
		tmp_h = geo_adapt(criteria__h,options);

		if ( (is_interactive == "true") || ( tmp_h.n_vertex() > 1.05*options.n_vertices_max))
		{
			// Show geometry, postprocess
			cout 	<< plotmtv << tmp_h;
			string feedback;
			cout	<< endl
					<< "*------------------------------------------------*" << endl
					<< "* Reread parameter file:                       r *" << endl
					<< "* Accept 					   y *" << endl
					<< "*------------------------------------------------*" << endl;
			cin >> feedback;
			//feedback ="y";
			if (feedback == "y")
			{
				flag = true;
			}
			else if (feedback == "r")
			{
			        //added by ali 12/16 Jun 2011
//				xmlparam->reparse_file();
			        load_xml_file();
				flag = false;
			}
		}
	}
	while (flag == false);

	omega_h = tmp_h;
	


    cout	<< endl
			<< "*------------------------------------------------*" << endl
			<< "* Adapted Mesh:                                  *" << endl
			<< "* hmin = " << omega_h.hmin() << endl
			<< "* #Elements = " << omega_h.size() << endl
			<< "*------------------------------------------------*" << endl;
    _adaptstep++;
	

	//ali 22 Oct 2011
	//save geo with full precision for restart
        omega_h.use_double_precision_in_saving();
	omega_h.save();

	_geoname = omega_h.name();
	cout 	<< "* - gemotry saved" << endl;

	// Reinitialize the FE geo, spaces, forms, fields
	multifield oldu(FEfields["u"]);
	multifield oldp(FEfields["p"]);
	multifield gamma;
	//added by ali 1 Jul 2011
	multifield oldT(FEfields["T"]);

	//added by ali 12/15 Jun 2011
	//xmlparam->get_childelement_string("SolverViscoplastic/type")
	std::string str_SolverViscoplastic_type;
	path_t solverVisc_type[] = {"SolverViscoplastic","type"};
	get_from_xml( solverVisc_type, &str_SolverViscoplastic_type );

	if (str_SolverViscoplastic_type=="Uzawa")
		gamma = FEfields["gamma"];

        FE_initialize();
	cout	<< "* - FE spaces reinitialised " << endl;

 	FEfields["u"].getfield() = interpolate(Vh, oldu.getfield() );
 	cout	<< "* - Reinterpolated u onto adapted mesh " << endl;
	FEfields["p"].getfield() = interpolate(Qh, oldp.getfield() );
	cout	<< "* - Reinterpolated p onto adapted mesh " << endl;
	if (str_SolverViscoplastic_type=="Uzawa")
	{
		FEfields["gamma"].getfield() = interpolate(Th, gamma.getfield() );
		cout	<< "* - Reinterpolated gamma onto adapted mesh " << endl;
		//added by ali 1 Jul 2011
                FEfields["T"].getfield() = interpolate(Th, oldT.getfield() );
                cout    << "* - Reinterpolated Tau onto adapted mesh " << endl;
	}		
	cout << "**************************************************" << endl << endl;

	//ali 25 Oct 2011
	//save for restart
	write_restart_files(0);

}

void solver_flowonset::set_basename()
{
    // Read basename from file:
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("basename") >> _basename;
  {
    path_t path[] = {"basename"};
    get_from_xml(path, &_basename);
  }

   _geoname = _basename;
   _dataname =_basename;

   cout << "\t Set basename to : " << _basename << endl;
  // Setup of the extensoin
  /* commented by ali 23 Oct 2011
    if (_timestep !=  0.)
    {
        char tmp[100];
        sprintf(tmp,"_rtime=%e",_timereal);
        _basename.append(tmp);
    }
    */
}



/**
 * @brief Set basename and write geometry, fields, and convergence rates (todo).
 *
 * @todo Write convergence rates.
 */
void solver_flowonset::write_all()
{
	cout 	<< "\t =========================================" << endl
			<< "\t = Write Files to disk:                   " << endl
			<< "\t -----------------------------------------" << endl;
    set_basename();
    write_geo();
    write_fields();
	cout 	<< "\t -----------------------------------------" << endl
			<< "\t - done" << endl
			<< "\t -----------------------------------------" << endl;
}


/**
 * @brief Write the fields to _basenema.fields.gz .
 *
 * Writes: u__h, p__h, gamma__h, gammadot__h, T__h
 */
void solver_flowonset::write_fields()
{
    // Open iorheostream
    orheostream out_u(_basename, "mfield");
	
	FEfieldsIterator iter;
	
	for (iter=FEfields.begin(); iter != FEfields.end();iter++)
	{
		out_u << "# geo= " << _basename <<".geo" << endl;
		out_u << "\n New field: \n\n";
    	out_u <<  rheo << catchmark(iter->first) << iter->second;	
	}	
    out_u.close();

}

/**
 * @brief Write the fields to _basenema.fields.gz .
 *
 * Writes: u__h, p__h, gamma__h, gammadot__h, T__h
 */
void solver_flowonset::write_fields(const string& name )
{
    // Open iorheostream
    orheostream out_u(name, "mfield");
	
	FEfieldsIterator iter;
	
	out_u << "# geo= " << _geoname <<".geo" << endl;
	out_u << "# dPdL= "<< *(monitors["dPdL"].end()-1);
        out_u << "# flowrate= Not computed";

	for (iter=FEfields.begin(); iter != FEfields.end();iter++)
	{
	    //added by ali 7 Jul 2011
	    if( iter->first=="f" || iter->first=="mf" )
	      continue;
	    out_u << "\n New field: \n\n";
    	    out_u <<  rheo << catchmark(iter->first) << iter->second;
	}	
    out_u.close();

    //added by ali 1 Aug 2011
    stringstream ss;
    ss << name << ".monitor";
    ofstream m(ss.str().c_str());
    assert( monitors["uL2"].size()==monitors["gL2"].size() );
    m << "deltaU   deltaGamma  (L2 norms)\n";
    for( size_t i=0; i<monitors["uL2"].size(); ++i ){
        m << monitors["uL2"][i] << '\t' << monitors["gL2"][i];
        m << '\n';
    }
    m << "\ndPdL iterations\n";
    for ( size_t i=0; i<monitors["dPdL"].size(); ++i ){
        m << monitors["dPdL"][i];
        m << '\n';
    }
    m.close();
    // ali 23 Oct 2011
    monitors["uL2"].clear();
    monitors["gL2"].clear();
}




/**
 * @brief Write the geometry file.
 */
void solver_flowonset::write_geo()
{
    orheostream out_geo(_basename, "geo");
    out_geo << omega_h;
    out_geo.close();
}

void solver_flowonset::write_geo(const string& name)
{
    orheostream out_geo(name, "geo");
    out_geo << omega_h;
    out_geo.close();
}



void solver_flowonset::solve()
{
	cout 	<< "\t==========================================\n"
    		<< "\t|Solver method                           |\n" 
    		<< "\t==========================================" << endl;
    int adapt_max=0;

    //added by ali 12/15 Jun 2011
//#    xmlparam->get_childelement_stream("gridadaption/gridadaptionsteps") >> adapt_max;
    path_t gridadap_steps_PT[] = {"gridadaption","gridadaptionsteps"};
    get_from_xml( gridadap_steps_PT, &adapt_max );

    cout	<< "\t| - Gridadaption steps = " << adapt_max << endl;

    assert( _adaptstep==0 );
    adapt_max -= _adaptstep;

    for (int i=0; i<= adapt_max; i++)
    {
        string algorithm;
        //added by ali 12/15 Jun 2011
        path_t solverVisc_type_PT[] = {"SolverViscoplastic","type"};
        get_from_xml( solverVisc_type_PT, &algorithm );

        std::string str_SolverViscoplastic_Model;
        path_t solverVisc_model_PT[] = {"SolverViscoplastic","Model"};
        get_from_xml( solverVisc_model_PT, &str_SolverViscoplastic_Model );

        if (algorithm == "Uzawa" & str_SolverViscoplastic_Model=="Bingham")
        {
            //added by ali 28 Jun 2011:
            ALG_Uzawa_test();
			stringstream tmp;
			tmp << _basename << "_gridadapt" << _adaptstep;
			_dataname = tmp.str();
			write_fields(_dataname);
            if( i!=adapt_max )
        	adapt_grid();
            else {
//	        viz1();
//	        plot_fp_fT();
//	        exec_shell("mkdir post && mv *.png *.gnu *.dat post/");
            }
        }
		else if (algorithm == "Uzawa" & str_SolverViscoplastic_Model!="Bingham")
		{
			VPREG_homotopy();
		}	
		else if (algorithm == "Regularised")
		{
			VPREG_homotopy();
		}
        //added by ali 7 Jul 2011
//        ++i_adapt;
    }
}


/*!
    \fn solver_flowonset::plot_fields()
	@brief Plot a menue and select the plottable vectorfields.
 */
void solver_flowonset::plot_fields(string title)
{
	bool whilebreak=true;	
	int choice;
	do
	{
		int i = 1;
		FEfieldsIterator fieldsiter;
		//map<string, multifield>::iterator fieldsiter;
		vector<string> selection;

		cout<< "\t|===========================================================|" << endl
			<< "\t| Graphical Output : " << title << endl
			<< "\t| MENU:  0 exit        100 options" << endl
			<< "\t|===========================================================|" << endl;
		for( fieldsiter = FEfields.begin(); fieldsiter != FEfields.end(); fieldsiter++ )
		{
			cout<< "\t| " << i << " " << fieldsiter->first << "\n";
			selection.push_back(fieldsiter->first);
			i++;
		}
		cout<< "\t|-----------------------------------------------------------|" << endl;
		cin >> choice;
		
		if (choice == 0)
			whilebreak = false;
		else if (0 < choice && choice <= i)
		{
			FEfields[selection[choice-1]].plot(selection[choice-1]);
		}
	}while(whilebreak);
}


/*!
    \fn solver_flowonset::analytic_solution()
	@brief Calculate the analytic solution of the problem. 

	@param amplitude	Amplitude of the wall perturbation
	@param Bn			Bingham number
	@param ywall		Coordinate of the wall
	@param dPdL			normalized pressure drop
	
	@todo Implement the asymptotic analysis
 */
void solver_flowonset::analytic_solution(Float amplitude, Float Bn, Float ywall,Float dPdL)
{
	if (amplitude == 0)
	{
		
		FEfields["analytic_u"] 		  = interpolate(Vh, analytic_straight_velocity(Bn,dPdL,1));
		FEfields["analytic_gammadot"] = interpolate(Th, analytic_straight_shearrate(Bn,dPdL,1));
		FEfields["analytic_u2shear"]  = FEforms["u2shear"]*FEfields["analytic_u"];
		FEfields["analytic_T"]		  = interpolate(Th, analytic_straight_stress(Bn,dPdL,1));
		multifield tt__h;
		FEfields["analytic_gamma"]	  = 0*FEfields["analytic_T"];
		tt__h 			 			  = FEfields["analytic_T"] + 1.0*FEfields["analytic_gammadot"];
		stress_to_gammadot(	tt__h,FEfields["analytic_gamma"],"Bingham",	Bn,2.0);
	}
	else
	{
		cout << " No analytic solution available " << endl;
	}
}


/*!
    \fn solver_flowonset::VPREG_homotopy()
	@brief Implements the homotopy loop for the flowrate problem. Also calls the grid adaption and interpolation.
	
	This algorithm does the following:
	- initial setup of the finite element fields and forms
	@todo Check for regularised problem: 
 */
void solver_flowonset::VPREG_homotopy()
{
    
	int r=0,rtot=0, maxregsteps;
	int solverstatus=0;
	Float Lambda, LambdaIn, LambdaTarget;
	string method;
	int adapt_max=0;

	//added by ali 12/15 Jun 2011
//#       xmlparam->get_childelement_stream("gridadaption/gridadaptionsteps") >> adapt_max;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/regparam_start") >> LambdaIn;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/regparam_end")   >> LambdaTarget;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/regparam_descend")>>method;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/maxregsteps")>>maxregsteps;
	{
	  path_t path[] = {"gridadaption","gridadaptionsteps"};
	  get_from_xml(path, &adapt_max);
	}
	{
	  path_t path[] = {"SolverViscoplastic","regparam_start"};
	  get_from_xml(path, &LambdaIn);
	}
	{
	  path_t path[] = {"SolverViscoplastic","regparam_end"};
	  get_from_xml(path, &LambdaTarget);
	}
	{
	  path_t path[] = {"SolverViscoplastic","regparam_descend"};
	  get_from_xml(path, &method);
	}
	{
	  path_t path[] = {"SolverViscoplastic","maxregsteps"};
	  get_from_xml(path, &maxregsteps);
	}
	
	cout 	<< "\t =========================================" << endl
			<< "\t = Homotopy method for Regularised method:" << endl
			<< "\t -----------------------------------------" << endl;
	do
	{
		Lambda = VPREG_Lambda(r, LambdaIn,LambdaTarget, method);
		_filenamecounters["Delta"]=Lambda;
		_filenamecounters["HomotopyStep"]=rtot;
		_modelparams["Delta"] = Lambda;
		_fluidmodel.set_params(_modelparams);
		
		cout	<< "\t =--------------------------------------------" << endl
				<< "\t = - Homotopy step: " << rtot << ", Lambda = " << Lambda << endl;


		//added by ali 12 Jun 2011
		//xmlparam->get_childelement_string("SolverViscoplastic/Model")
		std::string str_SolverViscoplastic_Model;
		path_t solverVisc_model_PT[] = {"SolverViscoplastic","Model"};
		get_from_xml( solverVisc_model_PT, &str_SolverViscoplastic_Model );

		if (str_SolverViscoplastic_Model!="Bingham")
		{
		        //added by ali 12/15 Jun 2011
		    //xmlparam->get_childelement_string("SolverViscoplastic/type")
		        std::string str_SolverViscoplastic_type;
		        path_t solverVisc_type_PT[] = {"SolverViscoplastic","type"};
		        get_from_xml( solverVisc_type_PT, &str_SolverViscoplastic_type );

			if (str_SolverViscoplastic_type=="Uzawa")
				ALG_Uzawa_test(Lambda);
			else
				solverstatus=Stokes_Piccard(false,Lambda);
		}
		else
			cout << "Incorrect model for the regularsed path" << endl;
				
		rtot++;r++;
		if (solverstatus==-1 && r>1)
		{
			r--;
			Float LambdaHat = VPREG_Lambda(r-1, LambdaIn,LambdaTarget, method);
			do
			{
				cout	<< "\t =--------------------------------------------" << endl
						<< "\t = - Homotopy Substep: " << rtot << ", Lambda = " << Lambda << endl;
				Lambda=(Lambda+LambdaHat)/2;
				_modelparams["Delta"] = Lambda;
				_fluidmodel.set_params(_modelparams);
				solverstatus=Stokes_Piccard(false, Lambda);
				rtot++;
			}while(solverstatus ==-1 && rtot<maxregsteps);	
		}
	
	// Declare Filenames and save file:
	stringstream tmp;
	tmp << _basename << "_Lambda=" << Lambda;
	_dataname = tmp.str();
	write_fields(_dataname);
	
	cout	<< "\t = - Wrote datafile: " << _dataname << endl;
	cout	<< "\t =--------------------------------------------" << endl;
	if (_adaptstep < adapt_max)
	{
		adapt_grid();
	}

	// plot_fields("Homotopy Loop");
		


	}while( rtot <maxregsteps && (Lambda>LambdaTarget || _adaptstep <= adapt_max) );
	
	cout 	<< "\t -----------------------------------------" << endl
			<< "\t = Homotopy sequence done" 
			<< "\t =========================================" << endl;

}



/**
 * @brief Calculates the regularisation param corresponding to homotopy step
 *
 * @param  r stepnumber
 * @param  LambdaIn Start value of reg. parameter
 * @param  LambdaTArget end value
 * @param  method decend method:
 * 			- "times0.1" 
 * @return Lambda
 */
Float  solver_flowonset::VPREG_Lambda(int r, Float LambdaIn,Float LambdaTarget,string method)
{
	Float Lambda;
    if (method == "times0.1" )
	{
		Lambda = pow(0.1,r)*LambdaIn;
	}
	if (Lambda < LambdaTarget)
		Lambda = LambdaTarget;
	return Lambda;
}



/**
 * @brief Piccard iteration for the flowrate problem: Setup of all necessary forms
 *
 * @param reg Regularisation parameter
 * @return 
	- 1 Solver converged inside to tolerance inside the maximal step limit
	- 2 Solver converges, but did not reach required tolerance inside the maximal step limit
	- -1 solver divergent
 */
int solver_flowonset::Stokes_Piccard(bool updateprecond, Float reg)
{
	cout 	<< "\t\t =============================================" << endl
			<< "\t\t Piccard Iteration for the nonlinear problem  " << endl
			<< "\t\t reg = " << reg <<endl
			<< "\t\t ---------------------------------------------" << endl;

	int miniter=0, maxiter,iter=0, StokesStatus=0;
	Float tol,toltarget=0.01;
	bool isflowrate=true;
	bool convergent=true;


	//added by ali 12/16 Jun 2011
	//# = xmlparam->get_childelement_string("SolverViscoplastic/Model");
	string model;
	{
	  path_t path[] = {"SolverViscoplastic","Model"};
	  get_from_xml(path, &model);
	}
	Float Bn;
//#	xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/nonlinmaxiter") >> maxiter;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/nonlinminiter") >> miniter;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/nonlintol") >> toltarget;
	{
	  path_t path[] = {"FluidProperties","Bn"};
	  get_from_xml(path, &Bn);
	}
	{
	  path_t path[] = {"SolverViscoplastic","nonlinmaxiter"};
	  get_from_xml(path, &maxiter);
	}
	{
	  path_t path[] = {"SolverViscoplastic","nonlinminiter"};
	  get_from_xml(path, &miniter);
	}
	{
	  path_t path[] = {"SolverViscoplastic","nonlintol"};
	  get_from_xml(path, &toltarget);
	}

	//added by ali 12/16 Jun 2011
	//xmlparam->get_childelement_string("subproblem")
	std::string str_subproblem;
	path_t subproblem_PT[] = {"subproblem"};
	get_from_xml( subproblem_PT, &str_subproblem );

	if (str_subproblem=="flowrate")
	{
		isflowrate=true;
		// Erase all unnecessary entries from the flowrate and the dPdL entry
		// in the monitor:
		if (monitors["dPdL"].size()>2)
			monitors["dPdL"].erase(monitors["dPdL"].begin(),monitors["dPdL"].end()-2);
		if (monitors["flowrate"].size()>2)
			monitors["flowrate"].erase(monitors["flowrate"].begin(),monitors["flowrate"].end()-2);
	} 
	else if (str_subproblem=="pressuredrop")
    	isflowrate=false;
	else
		cout << "Stokes_Piccard: unknown subproblem" << endl;
	
	monitors["PiccardResidual"].clear();

	// Test Step: Add analytic solutions for the straight channel
// 	analytic_solution(0, Bn, 1 , 10);
// 	FEfields["u"]=FEfields["analytic_u"];


	// Initialise the correct fields:
	FEfields["gammadot"] = FEforms["u2shear"]*FEfields["u"]; // Calculate rate of shear tensor
	FEfields["secinv_gammadot"] = secinv(FEfields["gammadot"]);
	secinvgammadot_to_viscosity(FEfields["viscosity"],FEfields["secinv_gammadot"], model, Bn, reg);
	//FEfields["viscosity"] = interpolate(FEfields["viscosity"].get_space(),_fluidmodel);
// 	plot_fields("Viscosity Calculation");
	do
	{
		cout << "\t\t  - Piccard Step " << iter << " of maximal " << maxiter << endl
			 << "\t\t    Store old u, viscosity and secinv_gammadot results for convergence tests\n";
		multifield oldu = FEfields["u"];	
		multifield oldviscosity = FEfields["viscosity"];
		multifield oldsecinv = FEfields["secinv_gammadot"];

		cout << "\t\t    Update weighted laplace form ..." ;
		FEforms["a"]=form(Vh,Vh,"2D_D",FEfields["viscosity"]);
		cout << "\t   ... done" << endl;
	
		//added by ali 12/16 Jun 2011
		//xmlparam->get_childelement_string("StokesSolver/type")
		std::string str_StokesSolver_type;
		path_t stokes_type_PT[] = {"StokesSolver","type"};

		if (get_string(stokes_type_PT)=="Uzawa")
		{
			Float augment;
			//added by ali 12/16 Jun 2011
//#			xmlparam->get_childelement_stream("StokesSolver/type") >> augment;
			path_t stokes_aug_PT[] = {"StokesSolver","type"};
			get_from_xml( stokes_aug_PT, &augment );

			FEforms["ar"] = FEforms["a"] + augment*trans(FEforms["b"])*FEforms["b"];
			if ( (updateprecond && iter % 50 == 0) || iter == 0 )
			{
				cout << "\t\t    Construct preconditioner ... ";
				fact = ldlt(FEforms["ar"].uu);
				// fact = 0*fact+EYE;
				//form tmp = form(Vh,Vh,"2D_D") + augment*trans(FEforms["b"])*FEforms["b"]; 
				//fact = ldlt(tmp.uu);
				updateprecond = false;
				cout << "\t   ... done" << endl;
			}
		}
		else
		{
			cout << "Incorrect solver type" << endl;
		}

		
		
		// Call solver
		if (str_subproblem=="flowrate")
		{
			StokesStatus = STOKES_flowrate();
		}
	
		// Update the fields depending on gammadot
		cout << "\t\t    Update gammadot, secinv, viscosity"  << endl;
		FEfields["gammadot"] = FEforms["u2shear"]*FEfields["u"]; // Calculate rate of shear tensor
		FEfields["secinv_gammadot"] = secinv(FEfields["gammadot"]);
		secinvgammadot_to_viscosity(FEfields["viscosity"],FEfields["secinv_gammadot"], model, Bn, reg);
		//FEfields["viscosity"] = interpolate(FEfields["viscosity"].get_space(),_fluidmodel);
		//plot_fields("Viscosity Calculation");

		Float tolu = norm(multifield(FEfields["u"]-oldu),"L2");
		Float tolv = norm(multifield(FEfields["viscosity"]-oldviscosity),"L2")/oldviscosity.max_abs();
		Float tolg = norm(multifield(FEfields["secinv_gammadot"]-oldsecinv),"L2");
		tol=max(tolu,tolv);

		cout << "\t\t    ||u-oldu||_L2           = " << tolu << endl
			 << "\t\t    ||visc-oldvisc||_L2/max = " << tolv << endl; 
		cout << "\t\t    ||secinv-oldsecinv||_L2 = " << tolg << endl; 
		
		monitors["PiccardResidual"].push_back(tol);
		iter++;	
	}while(iter < maxiter && (tol > toltarget || iter < miniter) && StokesStatus>0 && convergent);
// 	plot_fields("Piccard iterations end");
	
	cout 	<< "\t\t ---------------------------------------------  " << endl
			<< "\t\t Piccard Iterations done  " << endl
			<< "\t\t ---------------------------------------------" << endl;
	return 1;
}


void solver_flowonset::write_restart_files( const int niter )
{
  cout << "\t--------------------------------------------------------------------\n"
       << "\t|                        Saving for restart                         |"
       << endl;

  const size_t N = 10;
  if( N<monitors["dPdL"].size() )
        monitors["dPdL"].erase(monitors["dPdL"].begin(),monitors["dPdL"].end()-N);
  if( N<monitors["flowrate"].size() )
        monitors["flowrate"].erase(monitors["flowrate"].begin(),monitors["flowrate"].end()-N);

  std::ofstream o(restart_file);
  o << _adaptstep; o << '\n';
  o << niter; o << '\n';
  save_monitors(o);
  o.close();

  orheostream f("restart","mfield");
  f << setprecision( std::numeric_limits<Float>::digits10 );
  f << catchmark("u")     << FEfields["u"];
  f << catchmark("T")     << FEfields["T"];
  f << catchmark("gamma") << FEfields["gamma"];
  f.close();
  cout << "\t---------------------------------------------------------------------"
       << endl;
}


void solver_flowonset::save_monitors( std::ofstream& o )
{
  o << setprecision( std::numeric_limits<Float>::digits10 );
  Restart::save_std_vector( monitors["dPdL"], o, "dPdL" );
  Restart::save_std_vector( monitors["flowrate"], o, "flowrate" );
  Restart::save_std_vector( monitors["uL2"], o, "uL2" );
  Restart::save_std_vector( monitors["gL2"], o, "gL2" );
}


void solver_flowonset::load_monitors( std::ifstream& i )
{
  Restart::load_std_vector( monitors["dPdL"], i );
  Restart::load_std_vector( monitors["flowrate"], i );
  Restart::load_std_vector( monitors["uL2"], i );
  Restart::load_std_vector( monitors["gL2"], i );
}


//ali 21 Oct 2011
void solver_flowonset::restart( )
{
  std::ifstream info(restart_file);
  assert( info.is_open() );
  info >> _adaptstep;
  info >> uzawa_iter_restart;
  load_monitors(info);
  info.close();

  if( 0==_adaptstep )
    omega_h = geo(_basename+".geo");
  else {
    omega_h = geo(_basename+"-"+itos(_adaptstep)+".geo");
    omega_h.set_name(_basename);
    omega_h.data().set_serial_number(_adaptstep);
  }

  FE_initialize_spaces();
  FE_initialize_forms();
  FE_initialize_fields();

  irheostream i("restart","mfield");
  if(!i){
      cout << "Can not open restart file!---------------------\n";
      exit(0);
  }
  //note the order of reading should be the
  //same as what used when writing the file
  i >> catchmark("u")     >> FEfields["u"];
  i >> catchmark("T")     >> FEfields["T"];
  i >> catchmark("gamma") >> FEfields["gamma"];
  i.close();
}


void solver_flowonset::load_fields_for_postprocessing( const string& fieldfile )
{
 irheostream in;
 in.open(fieldfile,"mfield");
 in >> catchmark("u") >> FEfields["u"];
 in.close();

 in.open(fieldfile,"mfield");
 in >> catchmark("T") >> FEfields["T"];
 in.close();

 in.open(fieldfile,"mfield");
 in >> catchmark("gamma") >> FEfields["gamma"];
 in.close();

 in.open(fieldfile,"mfield");
 in >> catchmark("p") >> FEfields["p"];
 in.close();

 in.open(fieldfile,"mfield");
 in >> catchmark("gammadot") >> FEfields["gammadot"];
 in.close();

 in.open(fieldfile,"mfield");
 in >> catchmark("secinv_gamma") >> FEfields["secinv_gamma"];
 in.close();

 in.open(fieldfile,"mfield");
 in >> catchmark("secinv_gammadot") >> FEfields["secinv_gammadot"];
 in.close();
}


bool solver_flowonset::is_plug_broken( Float& wn, Float& ww )
{
  // assumptions:
  // 1- the shape of channel is symmetric one x and sym line is x=0
  // 2- min width of channel is on -L/2
  Float L;
  {
    path_t path[] = {"mesh","length"};
    get_from_xml(path,&L);
  }
  Myutil::set_Nan(&wn);
  Myutil::set_Nan(&ww);

  SegmentedLine line(point(-L/2,0),point(0,0),200,FEfields["p"].get_geo());

  enum {x,y};
  int in(0),iw(0);
  for( size_t i=0; i<line.size()-1; ++i )
  {
      rheolef::tensor gleft,grite;
      FEfields["gamma"].evaluate(line.point(i)  ,gleft);
      FEfields["gamma"].evaluate(line.point(i+1),grite);
      if( gleft(x,x)==0 && grite(x,x)!=0 ){
        wn = .5*L+line.point(i)[x];
        ++in;
      }
      if( gleft(x,x)!=0 && grite(x,x)==0 ){
        ww = -line.point(i+1)[x];
        ++iw;
      }
  }
  assert((in==0 && iw==0) || (in==1 && iw==1));
  return( in==1 );
}


struct segmentline_xiterator
{
  SegmentedLine const* l;
  size_t i;

  segmentline_xiterator( const SegmentedLine& line ): l(&line), i(0)
  {}

  enum {x,y};
  Float operator*() const
  {return l->point(i)[x];}

  void operator++()
  {++i;}
};

//interface adapter
struct Txx_iterator_adapter
{
  Field1D<rheolef::tensor> const* T;
  size_t i;

  Txx_iterator_adapter( const Field1D<rheolef::tensor>& t ): T(&t), i(0)
  {}

  enum {x,y};
  Float operator*()
  {return (T->at(i))(x,x);}

  void operator++()
  {++i;}
};

void solver_flowonset::plot_fp_fT()
{
  Float height, amplitude, L;
  {
    path_t path[] = {"mesh","height"};
    get_from_xml(path,&height);
  }
  {
    path_t path[] = {"mesh","amplitude"};
    get_from_xml(path,&amplitude);
  }
  {
    path_t path[] = {"mesh","length"};
    get_from_xml(path,&L);
  }

  Float wn,ww;
  Float lsample;
  if( is_plug_broken(wn,ww) ){
    lsample = 0.4*ww;
    printf("Broken plug with wn(%f) and ww(%f)\n",wn,ww);
  }
  else {
    lsample = L/6;
    printf("Intact plug\n");
  }

  enum {x,y};
  const Float max_width( height+amplitude );
  SegmentedLine symx( point(0,0),point(0,max_width),50,FEfields["p"].get_geo() );
  std::ofstream o("lsq.dat");
  o << "#y fP cP fT cT, Pressure and T show xfp(y) and xfT(y) curves\n";
  o << "#cP,cT are intercepts of least square fit\n";
  o << "#lsample " << lsample << '\n';
  Float mt_previous(0);
  Float y_dTxxdx_zero(0),dpdx_at_y_dTxxdx_zero;
  Myutil::set_Nan(&dpdx_at_y_dTxxdx_zero);
  for( size_t j=0; j<symx.size(); ++j )
  {
      const Float Y = symx.point(j)[y];
      o << Y <<' ';

      SegmentedLine sl(point(0,Y),point(-lsample,Y),20,FEfields["p"].get_geo());
      Postprocessing_FieldsOverLine R;
      put_fields_on_line(sl,R);

      segmentline_xiterator xiter(sl);
      Float mp,cp;
      Myutil::LinearLSQRFit lsq_fit;
      lsq_fit.set_npoints(sl.size());
      lsq_fit.leastsquare_linefit(xiter,R.P.begin(),&mp,&cp);
      o << mp <<' '<< cp <<'\t';

      Float mt,ct;
      Txx_iterator_adapter Titer(R.T);
      lsq_fit.leastsquare_linefit(xiter,Titer,&mt,&ct);
      o << mt <<' '<< ct <<'\n';

      if( mt*mt_previous<0 ){
        y_dTxxdx_zero = Y;
        dpdx_at_y_dTxxdx_zero = mp;
      }
      mt_previous = mt;
  }
  o.close();

  std::stringstream title;
  title << "'dTxxdx=0 @ y("<<y_dTxxdx_zero<<"), dpdx("<<dpdx_at_y_dTxxdx_zero<<")'\n";
  std::string basename("lsq");
  std::string fname("lsq.dat");
  std::string script =
  "set term png\n"
  "set output '"+basename+".png'\n"
  "set grid\n"
  "set key out\n"
  "set mxtics\n"
  "set mytics\n"
  "set title "+title.str()+
  "\n"
  "plot '"+fname+"' using 1:2 w l lw 2 t 'f(y) of P' \\\n"
  "    ,'"+fname+"' using 1:4 w l lw 2 t 'g(y) of T'\n";
  Postprocessing_FieldsOverLine::gen_gnuplot_script(basename,script);
  exec_shell("gnuplot "+basename+".gnu");
}


void solver_flowonset::put_fields_on_line( SegmentedLine& sline, Postprocessing_FieldsOverLine& r )
{
  using rheolef::Float;
  using rheolef::point;
  using rheolef::tensor;

  GradEvaluator GE( sline.spacing(), FEfields["p"].get_geo() );
  GE.remove_unvalid_points_of_line(sline);

  Field1D<MTinyVector<Float> >  gradp( sline, FEfields["p"], GE );
  Field1D<MTinyVector<point> >  gradv( sline, FEfields["u"], GE );
  Field1D<MTinyVector<tensor> > gradt( sline, FEfields["T"], GE );
  Field1D<MTinyVector<tensor> > gradg( sline, FEfields["gamma"], GE );

  r.GP = gradp;
  r.GV = gradv;
  r.GT = gradt;
  r.GGam = gradg;

  FieldEvaluator E;
  Field1D<Float> p( sline, FEfields["p"], E );
  Field1D<point> v( sline, FEfields["u"], E );
  Field1D<tensor> tau( sline, FEfields["T"], E );
  Field1D<tensor> gam( sline, FEfields["gamma"], E );

  r.P = p;
  r.V = v;
  r.T = tau;
  r.Gam = gam;


}


Postprocessing_FieldsOverLine solver_flowonset::horizn_output( const string& base, const string& title, const Float y )
{
  Float L;
  {
    path_t path[] = {"mesh","length"};
    get_from_xml(path,&L);
  }
  rheolef::point p1(-L/2,y);
  rheolef::point p2( L/2,y);
  SegmentedLine line(p1,p2,100,FEfields["p"].get_geo());

  Postprocessing_FieldsOverLine R;
  put_fields_on_line(line,R);
  write2file(base+".dat",line,R);
  Postprocessing_FieldsOverLine::viz_normalstress(base,title);
  exec_shell("gnuplot "+base+".gnu");
  return R;
}


void solver_flowonset::viz1( )
{
  double height, amplitude, L, Bn;
  {
    path_t path[] = {"mesh","height"};
    get_from_xml(path,&height);
  }
  {
    path_t path[] = {"mesh","amplitude"};
    get_from_xml(path,&amplitude);
  }
  {
    path_t path[] = {"mesh","length"};
    get_from_xml(path,&L);
  }
  {
    path_t path[] = {"FluidProperties","Bn"};
    get_from_xml(path,&Bn);
  }
  const double max_width( height+amplitude );


  enum {x,y};
  int N = 50;
  rheolef::point p1(0.,0.);
  rheolef::point p2(0.,max_width);
  string base;
  SegmentedLine line( p1, p2, N, FEfields["p"].get_geo() );
  Postprocessing_FieldsOverLine R;
  std::stringstream ss;
  ss << "deltah=" << 2*amplitude << ", Bn=" << Bn;

  base = "mid";
  put_fields_on_line(line,R);
  write2file(base+".dat",line,R);
  Postprocessing_FieldsOverLine::viz1(base,ss.str());
  exec_shell("gnuplot "+base+".gnu");

  //find height of yield surface and heigh of dead region (if exist)
  Float Hy(0.), Hd(0.);
  for( size_t i=0; i<line.size(); ++i )
  {
      if( (R.Gam(i))(x,y)!=0. )
      {
        Hy = (line.point(i))[y];
        //proceed to find H of dead region if exist
        for( int j=i+1; j<line.size(); ++j )
        {
            if( (R.Gam(j))(x,y)==0. ){
              Hd = (line.point(j))[y];
              break;
            }
        }
        break;
      }
  }
  assert(0.<Hy);
  const bool top_dead_region_exists = Hy<Hd;
  printf("Hy found: %f\n",Hy);
  if( top_dead_region_exists )
    printf("Hd found: %f\n",Hd);

  base = "Hy";
  std::stringstream ss2;
  ss2 << "Hy=" << std::setprecision(3) << Hy << ", " << ss.str();
  horizn_output("Hy",ss2.str(),Hy);

  base = "Hy2";
  ss2.str("");
  ss2 << "Hy/2=" << std::setprecision(3) << Hy/2 << ", " << ss.str();
  horizn_output("Hy2",ss2.str(),Hy/2);

  if( top_dead_region_exists )
    p2 = point(0,Hd);
  else
    p2 = point(0,max_width);
  const int M = 5;
  SegmentedLine l( point(0,Hy), p2, M, FEfields["p"].get_geo() );
  for( size_t i=1; i<M-1; ++i )
  {
      const double yc = (l.point(i))[y];
      std::stringstream t;
      t << "y=" << std::setprecision(3) << yc << ", " << ss.str();
      base = "y"+rheolef::itos(i);
      horizn_output(base,t.str(),yc);
  }


  if( top_dead_region_exists )
  {
    l.create_line( point(0,Hd), point(0,max_width), M );
    for( int i=0; i<M-1; ++i )
    {
      const double yd = (l.point(i))[y];
      base = "Hd"+itos(i);
      ss2.str("");
      ss2 << base+"=" << std::setprecision(4) << yd << ", " << ss.str();
      horizn_output(base,ss2.str(),yd);
    }
  }
}


template< class T >
void matlab( const field& f, const char* name, const geo& g )
{
  std::ofstream o(name);
  o << "# " << g.name() << '\n';
  for( geo::const_iterator_vertex vt=g.begin_vertex(); vt!=g.end_vertex(); ++vt )
  {
      T tmp;
      f.evaluate(*vt,tmp);
      GnuPlotWriter<T>::write(o,tmp);
      o << '\n';
  }
  o << '\n';
  o.close();

  string c = "gzip ";
  c.append(name);
  exec_shell(c);
}

// Feb 17, 2012
// writing fields for visualization in matlab
void solver_flowonset::dump_matlab( const char* fname )
{
  const geo& g = FEfields["p"].get_geo();
  exec_shell("mkdir matlab");

  matlab<Float>( FEfields["p"], "matlab/P", g );
  matlab<point>( FEfields["u"], "matlab/U", g );
  matlab<tensor>( FEfields["T"], "matlab/T", g );
  matlab<tensor>( FEfields["gamma"], "matlab/gamma", g );
  matlab<tensor>( FEfields["gammadot"], "matlab/gammadot", g );
  matlab<Float>( FEfields["secinv_gamma"], "matlab/secinv_gamma", g );
  matlab<Float>( FEfields["secinv_gammadot"], "matlab/secinv_gammadot", g );

  FEfields["secinv_T"] = secinv( FEfields["T"] );
  matlab<Float>( FEfields["secinv_T"], "matlab/secinv_T", g );
}
