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
 * \file solver_passiveparticle.cpp
 * \author Andreas Putz
 * \date 01-01-2008
 * \brief Class implementation of the steady state particle for viscoplastic (Bingham) fluids
 *
 * This program implements the main solver class for the wavychannel problem.
*/




#include "solver_passiveparticle.h"

//added by ali 12 Jun 2011
using namespace rheolef;

/**
 * @brief Constructor initialising the xmlparam
 * @param xml_parameters
 */
solver_passiveparticle::solver_passiveparticle( std::string const& fname_ )
:monitors()
{
    cout 	<< "\t==========================================\n"
    		<< "\t|Constructor Wavychennel                 |\n" 
			<< "\t==========================================" << endl;
    //added by ali 12 Jun 2011
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
	//added by ali 12 Jun 2011
//#	xmlparam->get_childelement_stream("ParticleProperties/velocityguess") >> pvel;
	monitors["ParticleVelocity"].push_back(pvel);
	
	FE_initialize();

	


	cout 	<< "\t|Constructor Wavychennel completed       |\n" 
			<< "\t==========================================" << endl << endl;
}


solver_passiveparticle::~solver_passiveparticle()
{}


void solver_passiveparticle::FE_initialize_geo()
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
int solver_passiveparticle::FE_initialize_spaces()
{
    // Initialize Spaces
	//cout << plotmtv << omega_h;
    string velspace,pspace,tspace,lmpspace;
    cout << "\t| - Setup spaces ..." << endl;

    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("SolverViscoplastic/FE_VelocitySpace") >> velspace;
//#    xmlparam->get_childelement_stream("SolverViscoplastic/FE_PressureSpace") >> pspace;
//#    xmlparam->get_childelement_stream("SolverViscoplastic/FE_StressSpace") >> tspace;
	
	cout << "\t|   VelSpace = " << velspace  << ", PSpace = " << pspace << ", StressSpace = " 	<< tspace << endl;

    Vh = space(omega_h,velspace,"vector");
    Qh = space(omega_h,pspace);
    Th = space(omega_h,tspace,"tensor");
	
    // Block for boundary conditions
    string geom;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("mesh/geom") >> geom;
    if (geom == "symxy")
    {
        Vh.block("top"); 		// No-Slip
	Vh.block("right");		// No-Slip
        Vh[0].block("bottom"); 	// Symmetry 	(u_x = 0)
        Vh[0].block("left");	// Symmetry 	(u_x = 0)
	Vh.block("particle");	// no horizontal component at particle surface
    }
    else if (geom == "full")
    {
        Vh.block("top"); 		// Block top
        Vh.block("bottom"); 	// Block bottom
        Vh.block("left");		// Inflow 		(uy = 0)
        Vh.block("right");		// Outflow 		(uy = 0)
	Vh.block("particle");	// no horizontal component at particle surface
    }
	
    cout << "\t|   Geometry type: " << geom << endl;
    cout << "\t|   *done\n" 
	 << "\t-------------------------------------------" << endl;
    return 1;
}


/*!
    \fn solver_passiveparticle::FE_initialize_fields()
 */
int solver_passiveparticle::FE_initialize_fields()
{
    cout << "\t| - Setup fields ..." << endl;
    string geom;
    string ViscoplSolver;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("SolverViscoplastic/type") >>  ViscoplSolver;

	if (ViscoplSolver == "Uzawa")
	{
		// Velocity contributions:
		FEfields["u"] 		 = multifield(Vh,0);
		FEfields["f"]		 = multifield(Vh,0);
		
		// Lagrange multiplier for the pressure and particle
		FEfields["p"] 		 = multifield(Qh,0);
		
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
		//added by ali 12 Jun 2011
//#		xmlparam->get_childelement_stream("FluidProperties/bodyforce") >>  magnitude;
		cout << "\t|    Bodyforce = " << magnitude << endl;
		FEfields["f"][1]=magnitude;
		FEfields["mf"] = FEforms["m"]*FEfields["f"];
	}
	
	
	// Initialise particle velocity:
	double velocity = *(monitors["ParticleVelocity"].end()-1);
	cout << "\t|    velocity = " << velocity << endl;
	FEfields["u"][1]["particle"]=velocity;
    return 1;
}


/*!
    \fn solver_passiveparticle::ALG_Uzawa()
	@brief Uzawa Algorithm for the augmented Lagrangian problem
	@todo Add description
 */
int solver_passiveparticle::ALG_Uzawa()
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
    string stokes_alg;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
//#    xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//#    xmlparam->get_childelement_stream("SolverViscoplastic/iterations") >> iter_uzawa_max;
//#    xmlparam->get_childelement_stream("StokesSolver/type") >> stokes_alg;
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
    \fn solver_passiveparticle::ALG_Uzawa_test()
	@brief Uzawa Algorithm for the augmented Lagrangian problem. Use the analytic solution.
	@todo Add description
 */
void solver_passiveparticle::ALG_Uzawa_test( Float Lambda )
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
    string stokes_alg;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
//#    xmlparam->get_childelement_stream("FluidProperties/Bn") >> Bn;
//#    xmlparam->get_childelement_stream("SolverViscoplastic/iterations") >> iter_uzawa_max;
//#    xmlparam->get_childelement_stream("StokesSolver/type") >> stokes_alg;
//#    xmlparam->get_childelement_stream("SolverViscoplastic/tolerance") >> tolerance;

    //added by ali 12 Jun 2011
    bool if_subproblem; // xmlparam->parse_childelement("subproblem")
    if (if_subproblem)
    {
        string subproblem;
//#        xmlparam->get_data_stream() >> subproblem;
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
		//added by ali 12 Jun 2011
		std::string str_SolverViscoplastic_Model; // xmlparam->get_childelement_string("SolverViscoplastic/Model")
		if (str_SolverViscoplastic_Model =="Bingham")
		{
			cout << "Bingham Model" << endl;
			stress_to_gammadot(	tt__h,
							FEfields["gamma"],"Bingham",
							Bn,1.0+r_Bi);
		}
		else if (str_SolverViscoplastic_Model =="Simple")
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
		// Convergence:
		FEfields["secinv_gammadot"]=secinv(FEfields["gammadot"]);
		FEfields["secinv_gamma"]=secinv(FEfields["gamma"]);
		Float tolu = norm(multifield(FEfields["u"]-oldu),"L2");
		Float toldgamma = norm(multifield(FEfields["secinv_gamma"]-FEfields["secinv_gammadot"]),"L2");
		tol=max(tolu,toldgamma);

		cout << "\t\t    ||u-oldu||_L2           = " << tolu << endl
			 << "\t\t    ||gamma-gammadot||_L2   = " << toldgamma << endl; 
		iter_uzawa++;
    }while(iter_uzawa<=iter_uzawa_max && tol > tolerance);
	
	//plot_fields("Field Output after all Uzawa steps");
	//string tt;
	//cin >> tt;
}

/*!
    \fn solver_passiveparticle::FE_initialize_forms()
	@brief Initialise the FE operator forms
 */
int solver_passiveparticle::FE_initialize_forms()
{
	cout << "\t| - Setup forms ..." << endl;
    /////////////////////////////////////////////////////////////////////////
    // Setup the integral forms (new way)
    /////////////////////////////////////////////////////////////////////////


	
	string type;
	//added by ali 12 Jun 2011
//#	xmlparam->get_childelement_stream("SolverViscoplastic/type") >> type;
    if (type == "Uzawa")
    {
		cout << "\t|   Viscoplasitic Solver: " << type;
		
		form a_h(Vh, Vh, "2D_D");
		form shearrate_tensor_h(Vh,Th,"2D");	// $int_\Omega \bs \shearr(\mb u):T dv$
		form inv_mass(Th,Th,"inv_mass"); 		// Inverse mass operator for 2 tensors
		form mass_h(Th,Th,"mass");				// Tensorial mass form for L2 norm
		
		
		// Mass forms for the rhs terms (velocity is test function)
		FEforms["m"] 		= form(Vh,Vh,"mass");		// Scalar product \int_Omega f \cdot v
		

		
		// Mass forms on the particle ( particle Lagrange multiplier is test function)
		
		// Constraint Matrix for the pressure
		FEforms["b"] = form(Vh, Qh, "div");
		FEforms["b"] = -FEforms["b"]; 			// grad, divergence operator
		
		FEforms["DIV_ThVh"] = -1. * trans(shearrate_tensor_h);
		FEforms["u2shear"]	= inv_mass*shearrate_tensor_h;

	//added by ali 12 Jun 2011
	bool if_SolverViscoplastic_augment; // xmlparam->parse_childelement("SolverViscoplastic/augment")
    	if (if_SolverViscoplastic_augment)
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
		        bool if_StokesSolver_type; //xmlparam->parse_childelement("StokesSolver/type")
			if (if_StokesSolver_type)
			{
			        //added by ali 12 Jun 2011
//#				xmlparam->get_data_stream() >> type;
				cout << "\t|   Stokes Solver: " << type;
				if (type == "Uzawa" || type == "urm_abtbc")
				{
				        //added by ali 12 Jun 2011
				        bool if_StokesSolver_augment; //xmlparam->parse_childelement("StokesSolver/augment")
					if (if_StokesSolver_augment)
					{
						double r_St;
						//added by ali 12 Jun 2011
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
int solver_passiveparticle::FE_initialize()
{
    cout 	<< "\t******************************************\n"
    		<< "\t|FE initialize                           |" << endl
    		<< "\t******************************************\n";
	FE_initialize_geo();
    FE_initialize_spaces();
    FE_initialize_forms();
	FE_initialize_fields();
	cout 	<< "\t|FE initialized                          |" << endl
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
int solver_passiveparticle::STOKES_Subproblem()
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
			//added by ali 12 Jun 2011
//#			xmlparam->get_childelement_stream("SolverViscoplastic/augment") >> r_Bi;
			tmpfield =  0.5*FEforms["DIV_ThVh"]*(FEfields["T"]-1*r_Bi*FEfields["gamma"]);
		}

	do
	{
		// Solve the Stokes problem:
		cout	<<"\t\t\t   - Stokes subproblem step " << iter +1
				<<" of " << itermax << endl;
		
		// Calculate the right hand side of the problem:
		field forcefield = FEfields["mf"];
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
int solver_passiveparticle::STOKES_Uzawa(field & rhs__h)
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
				FEfields["u"].u, FEfields["p"].u, rhs__h.u - (FEforms["ar"].ub*FEfields["u"].b),  	-(FEforms["b"].ub*FEfields["u"].b), 
				r_Stokes, iter_Stokes_max, tol_Stokes);
	
	// Reconstruct pressure and particle Lagrange multiplier:
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
field solver_passiveparticle::STOKES_rhs(field& mf__h, float r_Bi)
{
	//cout << "\t|\t Construct rhs";
    field mdiv__h = FEforms["DIV_ThVh"]*(FEfields["T"]-1*r_Bi*FEfields["gamma"]);
	//cout << plotmtv << mdiv__h;
    
	return 0.5*mdiv__h+FEfields["f"];
}



/**
 * @brief Adapt the grid
 *
 * This method uses either the shear rate tensor of the Lagrange multiplier
 * to conduct (bamg: anisotropic) grid refinement.
 * Also reinterpolates the velocity
 * @todo Reread XML file
 */
void solver_passiveparticle::adapt_grid()
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
			cout 	<< plotmtv << mayavi << tmp_h;
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
	//added by ali 12 Jun 2011
	std::string str_SolverViscoplastic_type; //xmlparam->get_childelement_string("SolverViscoplastic/type")
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

void solver_passiveparticle::set_basename()
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
void solver_passiveparticle::write_all()
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
void solver_passiveparticle::write_fields()
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
void solver_passiveparticle::write_fields(const string& name )
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

}




/**
 * @brief Write the geometry file.
 */
void solver_passiveparticle::write_geo()
{
    orheostream out_geo(_basename, "geo");
    out_geo << omega_h;
    out_geo.close();
}

void solver_passiveparticle::write_geo(const string& name)
{
    orheostream out_geo(name, "geo");
    out_geo << omega_h;
    out_geo.close();
}



void solver_passiveparticle::solve()
{
	cout 	<< "\t==========================================\n"
    		<< "\t|Solver method                           |\n" 
			<< "\t==========================================" << endl;
    int adapt_max=0;
    //added by ali 12 Jun 2011
//#    xmlparam->get_childelement_stream("gridadaption/gridadaptionsteps") >> adapt_max;
      cout << "\t| - Gridadaption steps = " << adapt_max << endl;
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
			//write_all();
        	adapt_grid();
        }
    }
}


/*!
    \fn solver_passiveparticle::plot_fields()
	@brief Plot a menue and select the plottable vectorfields.
 */
void solver_passiveparticle::plot_fields(string title)
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



