//
// C++ Implementation: make_grid
//
// Description: 
//
//
// Author: Andreas Putz <putza@math.ubc.ca>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "rheolef.h" 
#include "rheolef/uzawa_abtb.h"
#include <string>
#include <math.h>
using namespace std;


//added by ali 11 Jun 2011
using rheolef::Float;
using rheolef::ftos;
using rheolef::itos;



struct commandline{
  Float length;
  Float height;
  Float radius;
  Float amplitude;
  int subdivisions;
  string geom;
};

int main_a(int argc, char**argv) {
commandline cmd;
cmd.length = 20;
cmd.height  = 5;
cmd.radius = 1;
cmd.geom   = "symcos";
cmd.amplitude=-0.05;
cmd.subdivisions=100;
for (int i=1; i<argc; i++){
	if (!strcmp(argv[i],"-L") && ++i<argc)
		cmd.length=atof(argv[i]);
	else if (!strcmp(argv[i],"-H") && ++i<argc)
		cmd.height=atof(argv[i]);
	else if (!strcmp(argv[i],"-amp") && ++i<argc)
		cmd.amplitude=atof(argv[i]);
	else if (!strcmp(argv[i],"-geom") && ++i<argc)
		cmd.geom=argv[i];
	else if (!strcmp(argv[i],"-sub") && ++i<argc)
		cmd.subdivisions=atoi(argv[i]);
}
if (cmd.geom=="symcos"){
	string basename;
	basename = "symchannel_cos_L"+ftos(cmd.length);
	basename += "H"+ftos(abs(cmd.amplitude));
	// Create bamg file
	string bname = basename + "_g.msh";
	int subdivisions = 100;
	ofstream outf(bname.c_str());
	outf << "MeshVersionFormatted 0" << endl
		 << "AngleOfCornerBound 46" << endl
		 << "Dimension 2" << endl
		 << "Vertices " << itos(subdivisions+3)
		 << endl;
	outf << "-10 0 1" << endl;
	Float xi=-10;
	Float dx=cmd.length/subdivisions;
	Float pi=3.141592653589793238512808959406186204433;
	for(int i=0;i<=subdivisions;i++){
	  xi=-10 + i*dx;
	  outf  << xi
			<< " "
			<< 1+cmd.amplitude*cos((2*pi*xi)/cmd.length)
			<< " 1" 
			<< endl;
	}
	outf << "10 0 1" << endl;
	outf	<< "Edges "
			<< itos(subdivisions+3) 
			<< endl
			<< "1 2 1"
			<< endl;
	for(int i=2;i<=subdivisions+1;i++){
	  outf 	<< i << " "<< i+1 << " "
			<< "2" << endl;
	}
	outf << subdivisions+2 << " "
		 << subdivisions+3 <<" 3" << endl
		 << subdivisions+3 << " "
		 << "1 " <<" 4" << endl;
	outf << "hVertices" << endl;
	outf << "0.25 ";
	for(int i=1;i<=subdivisions;i++){
	  outf << cmd.length/subdivisions << " ";
	}
	outf << "0.25 " << "0.25 ";
	outf.close();
	string command ="bamg ";
	command += " -g " + bname;
	command += " -o " + basename + ".bamg";
	system (command.c_str());
	// Create domain list file:
	bname = basename + ".dmn";
	outf.open(bname.c_str());
	outf << "EdgeDomainNames" << endl
		 << "4" << endl
		 << "left" << endl
		 << "top" << endl
		 << "right" << endl
		 << "bottom";
	outf.close();
	
	// Transform into geo file
	command = "bamg2geo ";
	command += basename + ".bamg ";
	command += basename + ".dmn";
	command += " > " + basename + ".geo";
	system (command.c_str());

	// Show geo
	command = "geo " + basename + ".geo";
	system(command.c_str());
}
if (cmd.geom=="symsymcos"){
	string basename;
	basename = "symsymchannel_cos_L"+ftos(cmd.length);
	basename += "H"+ftos(abs(cmd.amplitude));
	// Create bamg file
	string bname = basename + "_g.msh";
	int subdivisions = 100;
	ofstream outf(bname.c_str());
	outf << "MeshVersionFormatted 0" << endl
		 << "AngleOfCornerBound 46" << endl
		 << "Dimension 2" << endl
		 << "Vertices " << itos(subdivisions+3)
		 << endl;
	outf << "-10 0 1" << endl;
	Float xi=-10;
	Float dx=cmd.length/subdivisions;
	Float pi=3.141592653589793238512808959406186204433;
	for(int i=0;i<=subdivisions;i++){
	  xi=-10 + i*dx;
	  outf  << xi
			<< " "
			<< 1-cmd.amplitude*cos((pi*xi)/cmd.length)
			<< " 1" 
			<< endl;
	}
	outf << "0 0 1" << endl;
	outf	<< "Edges "
			<< itos(subdivisions+3) 
			<< endl
			<< "1 2 1"
			<< endl;
	for(int i=2;i<=subdivisions+1;i++){
	  outf 	<< i << " "<< i+1 << " "
			<< "2" << endl;
	}
	outf << subdivisions+2 << " "
		 << subdivisions+3 <<" 3" << endl
		 << subdivisions+3 << " "
		 << "1 " <<" 4" << endl;
	outf << "hVertices" << endl;
	outf << "0.05 ";
	for(int i=1;i<=subdivisions;i++){
	  //outf << cmd.length/subdivisions << " ";
	  outf << "0.03" << " ";
	}
	outf << "0.03 " << "0.03 ";
	outf.close();
	string command ="bamg ";
	command += " -g " + bname;
	command += " -o " + basename + ".bamg";
	system (command.c_str());
	// Create domain list file:
	bname = basename + ".dmn";
	outf.open(bname.c_str());
	outf << "EdgeDomainNames" << endl
		 << "4" << endl
		 << "left" << endl
		 << "top" << endl
		 << "right" << endl
		 << "bottom";
	outf.close();
	
	// Transform into geo file
	command = "bamg2geo ";
	command += basename + ".bamg ";
	command += basename + ".dmn";
	command += " > " + basename + ".geo";
	system (command.c_str());

	// Show geo
	command = "geo " + basename + ".geo";
	system(command.c_str());
}
else if (cmd.geom=="symsymsphere"){
	string basename;
	basename = "symsymsphere_L"+ftos(cmd.length);
	basename += "_H"+ftos(abs(cmd.height));
	basename += "_R"+ftos(abs(cmd.radius));
	// Create bamg file
	string bname = basename + "_g.msh";
	int subdivisions = cmd.subdivisions;
	ofstream outf(bname.c_str());
	outf << "MeshVersionFormatted 0" << endl
		 << "AngleOfCornerBound 46" << endl
		 << "Dimension 2" << endl
		 << "Vertices " << itos(subdivisions+3)
		 << endl;
	Float pi=3.141592653589793238512808959406186204433;
	// Start with the circle:
	Float dphi=pi/(2*subdivisions);
	Float phi = pi/2;
	outf << "0 " << cmd.radius << " 5" << endl;
	for(int i=1;i< subdivisions-1;i++){
	  outf  << cmd.radius*cos(pi/2-i*dphi)
			<< " "
			<< cmd.radius*sin(pi/2-i*dphi)
			<< " 5" 
			<< endl;
	}
	outf << cmd.radius << " 0" << " 5" << endl;
	outf << cmd.length << " 0 2" << endl;
	outf << cmd.length << " " << cmd.height << " 3" << endl;
	outf << "0 " << cmd.height << " 4" << endl;
	
	// Write Edge information:
	outf	<< "Edges "
			<< itos(subdivisions+3) 
			<< endl;
	
	for(int i=1;i<subdivisions;i++){
	  outf 	<< i << " "<< i+1 << " "
			<< "5" << endl;
	}
	outf << subdivisions << " "
		 << subdivisions+1 <<" 1" << endl
		 << subdivisions+1 << " "
		 << subdivisions+2 << " 2" << endl
		 << subdivisions+2 << " "
		 << subdivisions+3 << " 3" << endl
		 << subdivisions+3 << " "
		 << 1              << " 4" << endl;
	
	outf << "hVertices" << endl;
	for(int i=1;i<=subdivisions+1;i++){
	  outf << 2*sin(dphi/2) << " ";
	}
	outf << "0.25 " << "0.25 ";
	outf.close();
	
	string command ="bamg ";
	command += " -g " + bname;
	command += " -o " + basename + ".bamg";
	system (command.c_str());
	// Create domain list file:
	bname = basename + ".dmn";
	outf.open(bname.c_str());
	outf << "EdgeDomainNames" << endl
		 << "5" << endl
		 << "bottom" 	<< endl
		 << "right" 	<< endl
		 << "top"  		<< endl
		 << "left" 		<< endl
		 << "sphere";
	outf.close();
	
	// Transform into geo file
	command = "bamg2geo ";
	command += basename + ".bamg ";
	command += basename + ".dmn";
	command += " > " + basename + ".geo";
	system (command.c_str());

	// Show geo
	command = "geo " + basename + ".geo";
	system(command.c_str());
}
}
