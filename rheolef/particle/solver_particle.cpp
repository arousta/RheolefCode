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
 * \file solver_particle.cpp
 * \author Andreas Putz
 * \date 01-01-2008
 * \brief Class implementation of the steady state particle for viscoplastic (Bingham) fluids
 *
 * This program implements the main solver class for the wavychannel problem.
*/

#include "solver_particle.h"

//added by ali 12 Jun 2011
using namespace rheolef;




/**
 * @brief Constructor initialising the xmlparam
 * @param xml_parameters
 */
solver_particle::solver_particle( std::string const& filename )
:monitors()
{
    cout 	<< "\t==========================================\n"
    		<< "\t|Constructor Wavychennel                 |\n" 
			<< "\t==========================================" << endl;
//    xmlparam = &xml_params;
    _timestep=0;
    _timereal=0.;
    _adaptstep=0;
	string velspace;

	set_basename();

	// Define some parameters:
	
	cout << "\t| - Setup monitors for the particle problem: ... \n";
	cout << "\t|   Initialise particle velocity\n";
	Float pvel;
//	xmlparam->get_childelement_stream("ParticleProperties/velocityguess") >> pvel;

	
	// Setup Monitors:
	monitors["ParticleVelocity"].push_back(pvel);
	monitors["a(u,u)"].push_back(0);
	monitors["j(u)"].push_back(0);
	monitors["L(u)"].push_back(0);
	
	FE_initialize();

	
	cout 	<< "\t|Constructor Wavychennel completed       |\n" 
		<< "\t==========================================" << endl << endl;
}


solver_particle::~solver_particle()
{}


void solver_particle::FE_initialize_geo()
{
	if (_adaptstep == 0)
	{
		cout << "\t| - Initialize GEO ... \n";
		//irheostream geostream(_basename,"geo");
		//omega_h = geo(_basename+".geo");
		omega_h = geo(_geoname+".geo");
		cout << "\t|    *done\n" 
			<< "\t-------------------------------------------" << endl;
	}
}


/**
 * @return flag
 * @brief Initialize the FE spaces (order, block DOF)
 */
int solver_particle::FE_initialize_spaces()
{
    // Initialize Spaces
	//cout << plotmtv << omega_h;
    string velspace,pspace,tspace,lmpspace;
    cout << "\t| - Setup spaces ..." << endl;

//    xmlparam->get_childelement_stream("SolverViscoplastic/FE_VelocitySpace") >> velspace;
//    xmlparam->get_childelement_stream("SolverViscoplastic/FE_PressureSpace") >> pspace;
//    xmlparam->get_childelement_stream("SolverViscoplastic/FE_StressSpace") >> tspace;
//    xmlparam->get_childelement_stream("ParticleProperties/FE_LagrangeMultSpace") >> lmpspace;
	cout << "\t|   VelSpace = " << velspace  << ", PSpace = " << pspace << ", StressSpace = " 	<< tspace << endl;
	cout << "\t|   Lagrange Multipler for particle motion = " << lmpspace << endl; 

	// Velocity Spaces
    Vh  = space(omega_h,velspace,"vector");
	Vh0 = space(omega_h,velspace,"scalar");
	Vh1 = space(omega_h,velspace,"scalar");
	// Pressure
    Qh = space(omega_h,pspace);
	// Stress / Gamma / Gammadot
    Th = space(omega_h,tspace,"tensor");
	// Lagrange multiplier for boundary
	LPh= space(omega_h,omega_h["particle"],lmpspace);//*space(omega_h,omega_h["particle"],lmpspace);


	// Block for boundary conditions
    string geom;
//    xmlparam->get_childelement_stream("mesh/geom") >> geom;


    if (geom == "symxy")
    {
		Vh.block("top");   Vh0.block("top");   Vh1.block("top");	// No-Slip
		Vh.block("right"); Vh0.block("right"); Vh1.block("right");	// No-Slip
		Vh[0].block("bottom");   Vh0.block("bottom"); 	// Symmetry 	(u_x = 0)
		Vh[0].block("left");     Vh0.block("left"); 	// Symmetry 	(u_x = 0)
		Vh[0].block("particle"); Vh0.block("particle");	// no horizontal component at particle surface
		//LPh[0].block("particle_from_ellipsoid");	// no horizontal component at particle surface 
    }
    else if (geom == "full")
    {
        Vh.block("top"); 		// Block top
        Vh.block("bottom"); 	// Block bottom
        Vh.block("left");		// Inflow 		(uy = 0)
        Vh.block("right");		// Outflow 		(uy = 0)
		Vh[0].block("particle");	// no horizontal component at particle surface
		//LPh[0].block("particle_from_ellipsoid");	// no horizontal component at particle surface 
    }
	
	Bh = Qh*LPh;
	
	// Diagnostic Output
	show_spacestats(Vh, "\t|   ", "Vh");
	show_spacestats(Vh0,"\t|   ", "Vh0");
	show_spacestats(Vh1,"\t|   ", "Vh1");
	
	show_spacestats(Qh, "\t|   ", "Qh");
	show_spacestats(LPh,"\t|   ", "LPhh");
	show_spacestats(Bh, "\t|   ", "Bh");
	
	
	cout << "\t|   Geometry type: " << geom << endl;
    cout << "\t|   *done\n" 
		 << "\t-------------------------------------------" << endl;
    return 1;
}


/*!
    \fn solver_particle::FE_initialize_fields()
 */
int solver_particle::FE_initialize_fields()
{
	cout << "\t| - Setup fields ..." << endl;
    string geom;



    string ViscoplSolver;
//    xmlparam->get_childelement_stream("SolverViscoplastic/type") >>  ViscoplSolver;



	if (ViscoplSolver == "Uzawa")
	{
		// Velocity contributions:
		FEfields["u"] 		 = multifield(Vh,0);
		FEfields["fp"]		 = multifield(Vh,0);
		FEfields["ff"]		 = multifield(Vh,0);
		
		// Lagrange multiplier for the pressure and particle
		FEfields["p"] 		 = multifield(Qh,0);
		FEfields["lambda"]	 = multifield(LPh,0);
		FEfields["U"]	 	 = multifield(LPh,0);
		FEfields["b"]		 = multifield(Bh,0);
		
		FEfields["gamma"] 	 = multifield(Th,0);
		FEfields["gammadot"] = multifield(Th,0);
		FEfields["T"] 		 = multifield(Th,0);
	}
	
	else
	{
		cout << "\t|    No valid method for the viscoplastic solver found" << endl;
	}

	// Impose force field on the fluid	
	{
		double magnitude;
//		xmlparam->get_childelement_stream("FluidProperties/bodyforce") >>  magnitude;


		FEfields["ff"][1]=magnitude;
		FEfields["mf"] = FEforms["m"]*FEfields["ff"];
	}
	// Impose force on the particle
	{
		double magnitude;
//		xmlparam->get_childelement_stream("ParticleProperties/coefficient") >>  magnitude;


		FEfields["fp"][1]["particle"]=magnitude/LPh.n_unknown();
		
	}
	//Construct total force contribution:
	FEfields["tf"] = FEfields["fp"]+FEfields["mf"];
	
	// Initialise particle velocity:
	FEfields["u"][1]["particle"]=*(monitors["ParticleVelocity"].end()-1);
    return 1;
}


/*!
    \fn solver_particle::ALG_Uzawa()
	@brief Uzawa Algorithm for the augmented Lagrangian problem
	@todo Add description
 */
int solver_particle::ALG_Uzawa()
{
	bool FLAG_DEBUG_GRAPH = true;
    // ===========================================================
    // Other local parameters
    // ===========================================================
    int iter_uzawa=0;

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

//    xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
//    xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//    xmlparam->get_childelement_stream("SolverViscoplastic/iterations") >> iter_uzawa_max;
//    xmlparam->get_childelement_stream("StokesSolver/type") >> stokes_alg;
    // -----------------------------------------------------------


    // ===========================================================
    // Main Uzawa Loop
    // ===========================================================

    for (iter_uzawa=0;iter_uzawa<=iter_uzawa_max;iter_uzawa++)
    {
        cout << "\t|   * timestep = " << _timestep 
	     << ", Adapt # = " << _adaptstep
	     << ", Iteration # = " << iter_uzawa << endl;

        // =======================================================
        // Step 1 (Stokes Problem)
        // =======================================================
		cout << "\t|     - Step 1 (Stokes Problem)" << endl;

        field mrhs__h =	 STOKES_rhs(mf__h,r_Bi); // Assemble rhs
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
		stress_to_gammadot(	tt__h,
							FEfields["gammadot"],"Bingham",
							Bn,1+r_Bi);
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

    }
	string tt;
	cin >> tt;
    return 1;
}


/*!
    \fn solver_particle::ALG_Uzawa_test()
	@brief Uzawa Algorithm for the augmented Lagrangian problem. Use the analytic solution.
	@todo Add description
 */
void solver_particle::ALG_Uzawa_test( Float Lambda )
{
	int FLAG_DEBUG_GRAPH = -1;
	cout <<"\t| - ALG_Uzawa" << endl;
    // ===========================================================
    // Other local parameters
    // ===========================================================
    int iter_uzawa=0;

    multifield tt__h; // temporary field
	
    // -----------------------------------------------------------

    // ==========================================================
    // Set local parameters from the XML File
    // ==========================================================
    Float r_Bi, Bn,tol,tolerance;
    int iter_uzawa_max;
    bool isflowrate = false;
    string stokes_alg;

//    xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
//    xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//    xmlparam->get_childelement_stream("SolverViscoplastic/iterations") >> iter_uzawa_max;
//    xmlparam->get_childelement_stream("StokesSolver/type") >> stokes_alg;
//    xmlparam->get_childelement_stream("SolverViscoplastic/tolerance") >> tolerance;

	Float UOld=0.,UNew=0.;


	
    if (true /*xmlparam->parse_childelement("subproblem")*/)
    {
        string subproblem;
//        xmlparam->get_data_stream() >> subproblem;
    }
    // -----------------------------------------------------------

	// ==================================================================
	// Define the Finite element matrices needed for one Uzawa iteration
	// ==================================================================
	
	// 1. Setup of the Stokes problem
	if ( stokes_alg == "urm_abtbc")	
		{
			// Pressure Stabilisation
			FEforms["c"] =  omega_h.hmin()*form(Qh, Qh , "grad_grad");
			//std::vector<int>::size_type nrow = c_h.nrow();

			cout << "\t|   Stokes Solver: " << stokes_alg << endl
				 << "\t| \t- Stabilisation term:\n" 
				 << "\t| \t  (c_h).nrow    = " << FEforms["c"].nrow()    << ", (c_h).col     = " << FEforms["c"].ncol() << endl
				 << "\t| \t  (c_h.uu).nrow = " << FEforms["c"].uu.nrow() << ", (c_h.uu).col  = " << FEforms["c"].uu.ncol() << endl
				 << "\t| \t  (c_h.ub).nrow = " << FEforms["c"].ub.nrow() << ", (c_h.ub).col  = " << FEforms["c"].ub.ncol() << endl
				 << "\t| \t  (c_h.bu).nrow = " << FEforms["c"].bu.nrow() << ", (c_h.bu).col  = " << FEforms["c"].bu.ncol() << endl
				 << "\t| \t  (c_h.bb).nrow = " << FEforms["c"].bb.nrow() << ", (c_h.bb).col  = " << FEforms["c"].bb.ncol() << endl;
			//cout << plotmtv << c_h;
		}	


	
	// plot_fields("Fields before Uzawa iteration");
	iter_uzawa = 0;
    do 
    {
        cout << "\t|   * timestep = " << _timestep 
			 << ", Adapt # = " << _adaptstep 
			 << ", Iteration # = " << iter_uzawa << endl;

		multifield oldu(FEfields["u"]);
        // =======================================================
        // Step 1 (Stokes Problem)
        // =======================================================
		cout << "\t|     - Step 1 (Stokes Problem)" << endl;	
		
		int StokesStatus = STOKES_Subproblem();
		
		FEfields["gammadot"] = FEforms["u2shear"]*FEfields["u"];
		
		if(((iter_uzawa % FLAG_DEBUG_GRAPH) == 0 && FLAG_DEBUG_GRAPH > 0))
			{
				plot_fields("Step 1: Velocity and derived rate of shear");
			}
		
 
        // =======================================================
        // Step 2 (Viscoplastic Problem)
        // =======================================================
		cout << "\t|     - Step 2 (Viscoplastic Problem): " ;
        // The use of the temprary stress is propably not optimal,
        // so it should be removed by implementing a conversion
        // method in the tensorfield class

		
		tt__h = FEfields["T"] + r_Bi*FEfields["gammadot"];

		//added by ali 11 Jun 2011
		std::string  sovler_viscoplastic_model;
//		sovler_viscoplastic_model = xmlparam->get_childelement_string("SolverViscoplastic/Model");
		if (sovler_viscoplastic_model =="Bingham")
		{
			cout << "Bingham Model" << endl;
			stress_to_gammadot(	tt__h,
							FEfields["gamma"],"Bingham",
							Bn,1.0+r_Bi);
		}
		else if (sovler_viscoplastic_model =="Simple")
		{
			cout << "Simple Model" << endl;
			stress_to_gammadot(	tt__h,
							FEfields["gamma"],"Simple",
							Bn,1.0+r_Bi,Lambda);
		}
		else
		{
			cerr << "No fluid model specified";
		}

		if((iter_uzawa % FLAG_DEBUG_GRAPH) == 0 && FLAG_DEBUG_GRAPH > 0)
			{
				plot_fields("Step 2: new gammadot");
			}
        

        // -------------------------------------------------------

        // =======================================================
        // Step 3 (Update of Lagrange Mult)
        // =======================================================
		// Again set analytic solution:
		//gammadot__h = agamma__h; T__h = aT__h;gamma__h=agamma__h;
		cout << "\t|     - Step 3 (Lagrange Mult)" << endl;
		FEfields["T"] = FEfields["T"] + 
			r_Bi*(FEfields["gammadot"]-FEfields["gamma"]);
        // -------------------------------------------------------

   
        // -------------------------------------------------------
		if( (iter_uzawa % FLAG_DEBUG_GRAPH) == 0 && FLAG_DEBUG_GRAPH > 0)
			{
				plot_fields("Field Output after all Uzawa substeps");
			}
			
			
		// =======================================================
        // Step 4 (Calculate particle velocity)
        // =======================================================
		cout << "\t|     - Step 4 (Calculate particle velocity)" << endl;
		
		FEfields["U"] = FEforms["VhToLPh"]*FEfields["u"];
		cout << "\t|       U(min) = " << FEfields["U"].min() << endl;
		cout << "\t|       U(max) = " << FEfields["U"].max() << endl;
		cout << "\t|       |U|    = " << norm(FEfields["U"]) << endl;
		cout << "\t|       U(max)-U(min) = " << FEfields["U"].max()- FEfields["U"].min() << endl;
		
		UNew = (FEfields["U"].max() + FEfields["U"].min())/2.;
			
		// =======================================================
        // Step 5 (Calculate tolerances)
        // =======================================================
		cout << "\t|     - Step 5 (Calculate tolerances)" << endl;
		
		FEfields["secinv_gammadot"]=secinv(FEfields["gammadot"]);
		FEfields["secinv_gamma"]=secinv(FEfields["gamma"]);
		Float tolu = norm(multifield(FEfields["u"]-oldu),"L2");
		Float toldgamma = norm(multifield(FEfields["secinv_gamma"]-FEfields["secinv_gammadot"]),"L2");
		Float tolU = abs(UNew-UOld);
		tol=max(tolu,toldgamma);tol = max(tol,tolU);

		cout << "\t\t    ||u-oldu||_L2           = " << tolu << endl
			 << "\t\t    ||gamma-gammadot||_L2   = " << toldgamma << endl
			 << "\t\t    |UNew-UOld|             = " << tolU << endl;
		
		UOld=UNew;
		iter_uzawa++;
    }while(iter_uzawa<=iter_uzawa_max && tol > tolerance);
	
	//plot_fields("Field Output after all Uzawa steps");
	//string tt;
	//cin >> tt;
}

/*!
    \fn solver_particle::FE_initialize_forms()
	@brief Initialise the FE operator forms
 */
int solver_particle::FE_initialize_forms()
{
	cout << "\t| - Setup forms ..." << endl;
    /////////////////////////////////////////////////////////////////////////
    // Setup the integral forms (new way)
    /////////////////////////////////////////////////////////////////////////


	
	string type;
//#	xmlparam->get_childelement_stream("SolverViscoplastic/type") >> type;


    if (type == "Uzawa")
    {
		cout << "\t|    - Viscoplasitic Solver: " << type << endl;
		
		form a_h(Vh, Vh, "2D_D");
		form shearrate_tensor_h(Vh,Th,"2D");	// $int_\Omega \bs \shearr(\mb u):T dv$
		form inv_mass(Th,Th,"inv_mass"); 		// Inverse mass operator for 2 tensors
		form mass_h(Th,Th,"mass");				// Tensorial mass form for L2 norm
		show_formstats(a_h,"\t|      ","2D_D");
		
		// Traces
		cout << "\t|    - Traces onto particle boundary: " << endl;
		FEforms["Vh0ToLPh"]  = trace(Vh0,LPh);		// Trace of global velocity space onto the particle
		FEforms["Vh1ToLPh"]  = trace(Vh1,LPh);		// Trace of global velocity space onto the particle
		{
			form_manip tmp;
			tmp << size(1,2)
				<<  FEforms["Vh0ToLPh"] << FEforms["Vh1ToLPh"];
			FEforms["VhToLPh"] = form(Vh,LPh);
			tmp >> FEforms["VhToLPh"];
		}
		show_formstats(FEforms["Vh0ToLPh"], "\t|      ", "Vh0ToLPh",true);
		show_formstats(FEforms["Vh1ToLPh"], "\t|      ", "Vh1ToLPh",true);
		show_formstats(FEforms["VhToLPh"],  "\t|      ", "VhToLPh", true);
		
		// Create modified trace operator (only transmitts V)
		FEforms["VhToLPhV"] = mk_modtrace("\t|      "); 	// Trace of global velocity space onto the particle with only one velocity
		show_formstats(FEforms["VhToLPhV"], "\t|      ", "VhToLphV", true );
		
		// Mass forms for the rhs terms (velocity is test function)
		FEforms["m"] 		= form(Vh,Vh,"mass");		// Scalar product \int_Omega f \cdot v
		
		// Mass forms on the particle ( particle Lagrange multiplier is test function)
		{
			Float factor=1;
			//added by ali 12 Jun 2011
			bool stokesSolver_augmentPart; //xmlparam->parse_childelement("StokesSolver/augmentPart")
			if (stokesSolver_augmentPart)
			{
//#				xmlparam->get_data_stream() >> factor;
			}
			FEforms["mp"] 		= factor*( form(LPh,LPh,"mass") + form(LPh,LPh,"d_ds") );	// Scalar product \int_{\Omega_P} f \cdot v
			cout << ", r_part = " << factor << ", ";
		}
		
		// Constraint Matrix for the pressure
		FEforms["b1"] = form(Vh, Qh, "div");
		FEforms["b1"] = -FEforms["b1"]; 			// grad, divergence operator
		
		// Constraint Matrix for the particle velocity (need modification with a projector):
		FEforms["b2"] = form(Vh,LPh);
		FEforms["b2"] = FEforms["mp"]*FEforms["VhToLPh"] - FEforms["mp"]*FEforms["VhToLPhV"];
// 		FEforms["b2"].uu = FEforms["mp"].uu*FEforms["VhToLPh"].uu - FEforms["mp"].uu*FEforms["VhToLPhV"].uu;
// 		FEforms["b2"].ub = FEforms["mp"].ub*FEforms["VhToLPh"].ub - FEforms["mp"].uu*FEforms["VhToLPhV"].ub;
				
		// Construct global constraint matrix "b":
		form_manip b_manip;
		b_manip	<< size(2,1)
				<< FEforms["b1"]
				<< FEforms["b2"];
		FEforms["b"] = form(Vh,Bh);
		b_manip >> FEforms["b"];
		
		
		
		// Forms for step 2/3:
		FEforms["DIV_ThVh"] = -1. * trans(shearrate_tensor_h);
		FEforms["u2shear"]	= inv_mass*shearrate_tensor_h;

	//added by ali 12 Jun 2011
	bool SolverViscoplastic_augment; //xmlparam->parse_childelement("SolverViscoplastic/augment")
    	if (SolverViscoplastic_augment)
        {
    		double r_Bi;
//#            xmlparam->get_data_stream() >> r_Bi;
            a_h = 1.0*r_Bi*a_h;
	    cout << ", r_Bi: " << r_Bi << endl;
        }

		// Setup fields for the Stokes solver
    	        //added by ali 12 Jun 2011
    	        bool if_StokesSolver; //xmlparam->parse_childelement("StokesSolver")
		if (if_StokesSolver)
		{
		        //added by ali 12 Jun 2011
		        bool if_StokesSolver_type; // xmlparam->parse_childelement("StokesSolver/type")
			if (if_StokesSolver_type)
			{
//#				xmlparam->get_data_stream() >> type;
				cout << "\t|   Stokes Solver: " << type;
				if (type == "Uzawa" || type == "urm_abtbc")
				{
				        //added by ali 12 Jun 2011
				        bool if_StokesSolver_augment; // xmlparam->parse_childelement("StokesSolver/augment")
					if (if_StokesSolver_augment)
					{
						double r_St;
//#						xmlparam->get_data_stream() >> r_St;
						FEforms["ar"] = a_h + r_St*trans(FEforms["b"])*FEforms["b"];
						fact = ldlt(FEforms["ar"].uu);
						cout << ", r_St: " << r_St << endl;
					}
				}
				else
				{
					cerr << "No valid Stokes solver specifier" << endl;
				}
			}
		}
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
int solver_particle::FE_initialize()
{
    cout << "\t******************************************\n"
    	 << "\t|FE initialize                           |" << endl
    	 << "\t******************************************\n";
    FE_initialize_geo();
    FE_initialize_spaces();
    FE_initialize_forms();
    FE_initialize_fields();
    cout << "\t|FE initialized                          |" << endl
    	 << "\t******************************************\n";
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
int solver_particle::STOKES_Subproblem()
{
	int iter=0, itermax=1;
	
	int solverstatus; 
	bool convergence=true;
	
	cout	<< "\t\t\t Stokes subproblem solver" << endl;

	field tmpfield;

	//added by ali 12 Jun 2011
	std::string str_SolverViscoplastic_type; //xmlparam->get_childelement_string("SolverViscoplastic/type")
	if (str_SolverViscoplastic_type=="Uzawa")
		{
			Float r_Bi;
//#			xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
			tmpfield =  0.5*FEforms["DIV_ThVh"]*(FEfields["T"]-1*r_Bi*FEfields["gamma"]);
		}

	do
	{
		// Solve the Stokes problem:
		cout	<<"\t\t\t   - Stokes subproble step " << iter +1
				<<" of " << itermax << endl;
		
		// Calculate the right hand side of the problem:
		field forcefield = FEfields["tf"];

		//added by ali 12 Jun 2011
		//xmlparam->get_childelement_string("SolverViscoplastic/type")
		if (str_SolverViscoplastic_type == "Uzawa" )
			forcefield = forcefield + tmpfield;

		//added by ali 12 Jun 2011
		std::string str_StokesSolver_type; //xmlparam->get_childelement_string("StokesSolver/type")
		if (str_StokesSolver_type == "Uzawa" )
			STOKES_Uzawa( forcefield );

		// Write the particle velocity to a monitor:
		// TODO: Extract the particle solution and write it to file
		{
// 			Float pvel = FEfields["u"][1]["boundary"].u[0];
// 			monitors["ParticleVelocity"].push_back(pvel);
		}
		
		// Setup of the secant method if itermax > 0
		// Check for double entries:
		iter++;
		//plot_fields("Flowfields inside flowrate calculation");
	}while(iter<itermax );

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
int solver_particle::STOKES_Uzawa(field & rhs__h)
{
    double r_Stokes;
    int iter_Stokes_max;
    double tol_Stokes;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("StokesSolver/augment")    >> r_Stokes;
//#    xmlparam->get_childelement_stream("StokesSolver/iterations") >> iter_Stokes_max;
//#    xmlparam->get_childelement_stream("StokesSolver/tolerance")  >> tol_Stokes;

// 	cout << "DOF u: unknowns:" << FEfields["u"].n_unknown()
// 		 << ", blocked: " << FEfields["u"].n_blocked() << endl;
// 	cout << "DOF ar: nrow:" << FEforms["ar"].nrow()
// 		 << "DOF ar: ncol:" << FEforms["ar"].ncol() << endl;
// 	cout << "DOF ar.ub: nrow:" << FEforms["ar"].ub.nrow()
// 		 << "DOF ar.ub: ncol:" << FEforms["ar"].ub.ncol() << endl;

    uzawa_abtb (FEforms["ar"].uu, fact, FEforms["b"].uu,
				FEfields["u"].u, FEfields["b"].u, rhs__h.u - (FEforms["ar"].ub*FEfields["u"].b),  	-(FEforms["b"].ub*FEfields["u"].b), 
				r_Stokes, iter_Stokes_max, tol_Stokes);
	
	// Reconstruct pressure and particle Lagrange multiplier:
	FEfields["p"]=FEfields["b"][0];
	FEfields["lambda"][0]=FEfields["b"][1];
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
field solver_particle::STOKES_rhs(field& mf__h, float r_Bi)
{
	//cout << "\t|\t Construct rhs";
    field mdiv__h = FEforms["DIV_ThVh"]*(FEfields["T"]-1*r_Bi*FEfields["gamma"]);
	//cout << plotmtv << mdiv__h;
    
	return 0.5*mdiv__h+FEfields["tf"];
}



/**
 * @brief Adapt the grid
 *
 * This method uses either the shear rate tensor of the Lagrange multiplier
 * to conduct (bamg: anisotropic) grid refinement.
 * Also reinterpolates the velocity
 * @todo Reread XML file
 */
void solver_particle::adapt_grid()
{
    string is_interactive;
    string criteria;
    Float Bn;
    field criteria__h;
	geo tmp_h;


	cout 	<< "\n**************************************************" << endl
        	<< "* Adaptive step                                  *" << endl;
	bool flag = false;
	do
	{
		flag = true;
		//added by ali 12 Jun 2011
//#		xmlparam->reparse_file();
//#		xmlparam->get_childelement_stream("gridadaption/interactive") >> is_interactive;
//#    	        xmlparam->get_childelement_stream("gridadaption/criteria") >> criteria;
		adapt_option_type options;
//#		xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//#		xmlparam->get_childelement_stream("gridadaption/max_vertices") >> options.n_vertices_max;
//#		xmlparam->get_childelement_stream("gridadaption/bamg_options/hcoef") >> options.hcoef;
//#		xmlparam->get_childelement_stream("gridadaption/bamg_options/hmin") >> options.hmin;
//#		xmlparam->get_childelement_stream("gridadaption/bamg_options/hmax") >> options.hmax;
//#		xmlparam->get_childelement_stream("gridadaption/bamg_options/hratio") >> options.ratio;

		

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
			cout 	<< mayavi << tmp_h;
			string feedback;
			cout	<< endl
					<< "*------------------------------------------------*" << endl
					<< "* Reread parameter file:                       r *" << endl
					<< "* Accept 					                   y *" << endl
					<< "*------------------------------------------------*" << endl;
			cin >> feedback;
			//feedback ="y";
			if (feedback == "y")
			{
				flag = true;
			}
			else if (feedback == "r")
			{
			        //added by ali 12 Jun 2011
//#				xmlparam->reparse_file();
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
	

	omega_h.save();
	_geoname = omega_h.name();
	cout 	<< "* - gemotry saved" << endl;

	// Reinitialize the FE geo, spaces, forms, fields
	multifield oldu(FEfields["u"]);
	multifield oldp(FEfields["p"]);
	multifield gamma;

	//added by ali 12 Jun  2011
	std::string str_SolverViscoplastic_type; // xmlparam->get_childelement_string("SolverViscoplastic/type")
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
	}		
	cout << "**************************************************" << endl << endl;
}

void solver_particle::set_basename()
{
    // Read basename from file:
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("basename") >> _basename;
    // Setup of the extensoin
    if (_timestep !=  0.)
    {
        char tmp[100];
        sprintf(tmp,"_rtime=%e",_timereal);
        _basename.append(tmp);
    }
	_geoname = _basename;
	_dataname =_basename;
	
    cout << "\t Set basename to : " << _basename << endl;
}



/**
 * @brief Set basename and write geometry, fields, and convergence rates (todo).
 *
 * @todo Write convergence rates.
 */
void solver_particle::write_all()
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
void solver_particle::write_fields()
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
void solver_particle::write_fields(const string& name )
{
    // Open iorheostream
    orheostream out_u(name, "mfield");
	
	FEfieldsIterator iter;
	
	out_u << "# geo= " << _geoname <<".geo" << endl;
	//out_u << "# flowrate= "<< *(monitors["flowrate"].end()-1);

	for (iter=FEfields.begin(); iter != FEfields.end();iter++)
	{
		out_u << "\n New field: \n\n";
    	out_u <<  rheo << catchmark(iter->first) << iter->second;	
	}	
    out_u.close();
	
	for (iter=FEfields.begin(); iter != FEfields.end();iter++)
	{
		// Write corresponding file:
		out_u << "\n New field: \n\n" <<  iter->first;
		(iter->second).save( name + iter->first);
	}	
}



/**
 * @brief Write the stats
 *
 * Writes: U, a(u,u), j(u), L(u)
 */
void solver_particle::write_stats(const string& name )
{
    // Open iorheostream
	cout 	<< "\t =========================================" << endl
			<< "\t = Write data stats:                      " << endl
			<< "\t -----------------------------------------" << endl;
	
	orheostream out_u(name, "dat");
	
	FEfieldsIterator iter;
	{
		
		out_u << "U = " << ( FEfields["U"].max() + FEfields["U"].min() )/2 << endl;
	}
	// Write a(u,u);
	{
		form a_h(Vh, Vh, "2D_D");
		out_u << "a(u,u) = " << dot(FEfields["u"],a_h*FEfields["u"])<< endl;
	}
	// Write j(u);
	{
		space tmpspace;
		tmpspace = Th[0];
		form mass_h(tmpspace,tmpspace,"mass");
		field tmpfield(FEfields["secinv_gamma"].get_space(),1);
		out_u << "j(u)   = " << dot( tmpfield , mass_h*FEfields["secinv_gamma"] ) << endl;
	}
	// Write L(u);
	{
		out_u << "L(u)   = " << dot( FEfields["u"] , FEfields["tf"] ) << endl;
	}
	out_u.close();
	cout 	<< "\t -----------------------------------------" << endl
			<< "\t - done" << endl
			<< "\t -----------------------------------------" << endl;
}



/**
 * @brief Write the geometry file.
 */
void solver_particle::write_geo()
{
    orheostream out_geo(_basename, "geo");
    out_geo << omega_h;
    out_geo.close();
}

void solver_particle::write_geo(const string& name)
{
    orheostream out_geo(name, "geo");
    out_geo << omega_h;
    out_geo.close();
}



void solver_particle::solve()
{
	cout << "\t==========================================\n"
    	     << "\t|Solver method                           |\n"
	     << "\t==========================================" << endl;
    int adapt_max=0;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("gridadaption/gridadaptionsteps") >> adapt_max;
	cout	<< "\t| - Gridadaption steps = " << adapt_max << endl; 
    for (int i=0; i<= adapt_max; i++)
    {
        string algorithm;
        //added by ali 12 Jun 2011
//#        xmlparam->get_childelement_stream("SolverViscoplastic/type") >> algorithm;
        std::string str_SolverViscoplastic_Model; //xmlparam->get_childelement_string("SolverViscoplastic/Model")
        if (algorithm == "Uzawa" & str_SolverViscoplastic_Model=="Bingham")
        {
            ALG_Uzawa_test();
			stringstream tmp;
			tmp << _basename << "_gridadapt=" << _adaptstep;
			_dataname = tmp.str();
			write_fields(_dataname);
			write_stats(_dataname);
			//write_all();
        	adapt_grid();
        }
    }
}


/*!
    \fn solver_particle::plot_fields()
	@brief Plot a menue and select the plottable vectorfields.
 */
void solver_particle::plot_fields(string title)
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





/**
 * Create the modified trace operator onto the boundary
 * @return 
 */
form solver_particle::mk_modtrace( string ident )
{
	cout 	<< endl;
	cout 	<< ident << " Create modified trace onto boundary: " << endl;
			
	form tmpform(Vh,LPh);
	form t1(Vh1,LPh);
	
	{
		cout 	<< ident << "   - Create Vh1 part: UU " << endl;
		int ncol = FEforms["Vh1ToLPh"].uu.ncol();
		int nrow = FEforms["Vh1ToLPh"].uu.nrow();
		int nnz  = FEforms["Vh1ToLPh"].uu.nnz();
		
		cout	<< ident << "     Dimensions of Vh1 part:" << endl
				<< ident << "      - nrow =  " << nrow << endl
				<< ident << "      - ncol =  " << ncol << endl
				<< ident << "      - nnz =  " <<  nnz  << endl;
		
		csr<Float> tmpmatrix = FEforms["Vh1ToLPh"].uu;
		
		dns<Float> dense_ones(nrow,1);
		dense_ones=1;
		csr<Float> matrix_ones(dense_ones);
		
		csr<Float> matrix_zeros(nrow,nrow-1,0);
		matrix_zeros=0;
		
		csr<Float> matrix_map= hcat(matrix_ones,matrix_zeros);
		
		t1.uu = matrix_map*FEforms["Vh1ToLPh"].uu;
	}
	
	{
		cout 	<< ident << "   - Create Vh1 part: UB " << endl;
		int ncol = FEforms["Vh1ToLPh"].ub.ncol();
		int nrow = FEforms["Vh1ToLPh"].ub.nrow();
		int nnz  = FEforms["Vh1ToLPh"].ub.nnz();
		
		cout	<< ident << "     Dimensions of Vh1 part:" << endl
				<< ident << "      - nrow =  " << nrow << endl
				<< ident << "      - ncol =  " << ncol << endl
				<< ident << "      - nnz =  " <<  nnz  << endl;
		
		csr<Float> tmpmatrix = FEforms["Vh1ToLPh"].ub;
		
		csr<Float> matrix_zeros(nrow,ncol,0);
		matrix_zeros=0;
		
		t1.ub = matrix_zeros;
	}
	
	{
		cout 	<< ident << "   - Create Vh1 part: BU " << endl;
		int ncol = FEforms["Vh1ToLPh"].bu.ncol();
		int nrow = FEforms["Vh1ToLPh"].bu.nrow();
		int nnz  = FEforms["Vh1ToLPh"].bu.nnz();
		
		cout	<< ident << "     Dimensions of Vh1 part:" << endl
				<< ident << "      - nrow =  " << nrow << endl
				<< ident << "      - ncol =  " << ncol << endl
				<< ident << "      - nnz =  " <<  nnz  << endl;
		
		csr<Float> tmpmatrix = FEforms["Vh1ToLPh"].bu;

		
		csr<Float> matrix_zeros(nrow,ncol,0);
		matrix_zeros=0;
		
		
		t1.bu = matrix_zeros;
	}
	
	{
		cout 	<< ident << "   - Create Vh1 part: BB " << endl;
		int ncol = FEforms["Vh1ToLPh"].bb.ncol();
		int nrow = FEforms["Vh1ToLPh"].bb.nrow();
		int nnz  = FEforms["Vh1ToLPh"].bb.nnz();
		
		cout	<< ident << "     Dimensions of Vh1 part:" << endl
				<< ident << "      - nrow =  " << nrow << endl
				<< ident << "      - ncol =  " << ncol << endl
				<< ident << "      - nnz =  " <<  nnz  << endl;
		
		csr<Float> tmpmatrix = FEforms["Vh1ToLPh"].bb;
		
		
		csr<Float> matrix_zeros(nrow,ncol,0);
		matrix_zeros=0;
		
		t1.bb = matrix_zeros;
	}
	show_formstats(t1,ident+"   ","testform",true);
	
	form_manip tmp_manip;
	tmp_manip 	<< size(1,2)
				<< FEforms["Vh0ToLPh"] << t1;
	tmp_manip >> tmpform;
	return tmpform;
}


void solver_particle::show_formstats(const form & testform, string ident, string name, bool plots)
{
	cout 	<< ident << " Name of the form: " << name << endl
			<< ident << "  - nrow        = " << testform.nrow() << endl
			<< ident << "  - ncol        = " << testform.ncol() << endl;
	
	cout 	<< ident << "  Matrix: UU " << endl
			<< ident << "   - nrow " << testform.uu.nrow() << endl
			<< ident << "   - ncol " << testform.uu.ncol() << endl
			<< ident << "   - nnz  " << testform.uu.nnz() << endl;
	
	cout 	<< ident << "  Matrix: UB " << endl
			<< ident << "   - nrow " << testform.ub.nrow() << endl
			<< ident << "   - ncol " << testform.ub.ncol() << endl
			<< ident << "   - nnz  " << testform.ub.nnz() << endl;
	
	cout 	<< ident << "  Matrix: BU " << endl
			<< ident << "   - nrow " << testform.bu.nrow() << endl
			<< ident << "   - ncol " << testform.bu.ncol() << endl
			<< ident << "   - nnz  " << testform.bu.nnz() << endl;
	
	cout 	<< ident << "  Matrix: BB " << endl
			<< ident << "   - nrow " << testform.bb.nrow() << endl
			<< ident << "   - ncol " << testform.bb.ncol() << endl
			<< ident << "   - nnz  " << testform.bb.nnz() << endl;
	
	if (plots)
	{
		{
				// Write matrix to file:
			string bname = name + "_UU";
			orheostream out_u(bname, "mtx");
			out_u << matrix_market << testform.uu;
			out_u.close();
		}
		
		{
				// Write matrix to file:
			string bname = name + "_UB";
			orheostream out_u(bname, "mtx");
			out_u << matrix_market << testform.ub;
			out_u.close();
		}
		
		{
				// Write matrix to file:
			string bname = name + "_BU";
			orheostream out_u(bname, "mtx");
			out_u << matrix_market << testform.bu;
			out_u.close();
		}
		
		{
				// Write matrix to file:
			string bname = name + "_BB";
			orheostream out_u(bname, "mtx");
			out_u << matrix_market << testform.bb;
			out_u.close();
		}
	}
}

void solver_particle::show_spacestats(const space & testspace, string ident="", string spacename="")
{
	cout 	<< ident << " Name of the Space: " << spacename << endl
			<< ident << "  - size        = " << testspace.size() << endl
			<< ident << "  - n_blocked   = " << testspace.n_blocked() << endl
			<< ident << "  - n_unknown   = " << testspace.n_unknown() << endl
			<< ident << "  - dimension   = " << testspace.dimension() << endl
			<< ident << "  - n_components= " << testspace.n_component () << endl;
	int ncomp = testspace.n_component ();
	if(ncomp >1)
	{
		for(int i = 0; i < ncomp; i++)
		{
			cout 	<< ident << "   Component " <<  i << endl
					<< ident << "  - size        = " << testspace.size_component(i) << endl
					<< ident << "  - n_blocked   = " << testspace.n_blocked_component(i) << endl
					<< ident << "  - n_unknown   = " << testspace.n_unknown_component(i) << endl;
		}
	}
}
