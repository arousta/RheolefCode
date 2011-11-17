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
#include "mk_mesh.h"

#include "SimpleXML.h"
#include "addedUtil.h"
//ali Nov 2011
#include "FieldVisualization.h"


using namespace std;

//added by ali 11 Jun 2011
using rheolef::itos;

//added by ali 29 Jun 2011
mk_mesh::mk_mesh( char const* fname ):
xmldoc( new TiXmlDocument )
{
  if( !(xmldoc->LoadFile(fname)) ){
    cout << "Can not load the file " << fname << endl;
    exit(1);
  }

  params = xmldoc->RootElement()->FirstChildElement("parameters");
}


//added by ali 29 Jun 2011
mk_mesh::~mk_mesh()
{
  delete xmldoc;
}

/**
 * 
 * @brief Creates the meshed
 * @param xmlparam XML input file
 * @return flag
 */
int mk_mesh::bamg_makemesh( )
{
	string meshtype;
//	xmlparam.get_childelement_stream("problem") >> meshtype;
	path_t path[] = {"problem"};
	get_from_xml(path, &meshtype);

	if (meshtype=="wavychannel")
		bamg_wavychannel();
	else if (meshtype=="ellipsoid")
		bamg_ellipsoid();
	else if(meshtype=="squarechannel")
	        bamg_squarechannel();
	else if(meshtype=="trichannel")
	        bamg_trichannel();
	else
	{
		cout 	<< "Wrong problem supplied:" << endl
				<< "  You supplied: " << meshtype << endl
				<< "  Allowed values are: " << endl
				<< "     - wavychannel" << endl
				<< "     - ellipsoid" << endl
				<< "     - squarechannel" << endl
				<< "     - trichannel"  << endl;
	}
	return 1;
}


void mk_mesh::bamg_trichannel()
{
      string basename,filename,path;
      string geom;
      double length,height,amplitude,width,hmax;
      int subdivisions_top, nvertices;

      Float xi=0.,dx=0.;

      {
        path_t ppath[] = {"workdir"};
        get_from_xml(ppath, &path);
      }
      {
        path_t path[] = {"basename"};
        get_from_xml(path, &basename);
      }
      filename = path + "/" + basename;

      {
        path_t path[] = {"mesh","geom"};
        get_from_xml(path, &geom);
      }
      {
        path_t path[] = {"mesh","length"};
        get_from_xml(path, &length);
      }
      {
        path_t path[] = {"mesh","height"};
        get_from_xml(path, &height);
      }
      {
        path_t path[] = {"mesh","width"};
        get_from_xml(path,&width);
      }
      {
        path_t path[] = {"mesh","amplitude"};
        get_from_xml(path, &amplitude);
      }
      {
        path_t path[] = {"mesh","subdivisions_top"};
        get_from_xml(path, &subdivisions_top);
      }
      {
        path_t path[] = {"mesh","hmax"};
        get_from_xml(path, &hmax);
      }

      cout << "\n\t\tWorkdir              = " << path << endl
              << "\n\t\tbasename              = " << basename << endl
              << "\n\t\tmesh/geom             = " << geom << endl
              << "\n\t\tmesh/length           = " << length << endl
              << "\n\t\tmesh/height           = " << height << endl
              << "\n\t\tmesh/amplitude        = " << amplitude << endl
              << "\n\t\tmesh/hmax             = " << hmax << endl
              << "\n\t\tmesh/subdivsions_top  = " << subdivisions_top << endl;
      cout << "...done\n" << endl;

      if ("symy" == geom)
      {
              dx=length/subdivisions_top; //stepsize
              nvertices = 7;
      }
      else if ("symxy" == geom)
      {
              length = length*2;
              dx=length/(2.*subdivisions_top);
              nvertices = 5;
      }
      else
      {
              cerr << "Wrong geometry type !!!" << endl;
              std::exit(-1);
      }

      string bname = filename + "_g.msh";
      ofstream outf(bname.c_str());
      outf << "MeshVersionFormatted 0\n"
           << "AngleOfCornerBound 46\n"
           << "Dimension 2\n"
           << "Vertices " << nvertices
           << '\n';

      double ych1 = height-amplitude;
      double ych2 = height+amplitude;
      xi= -length/2.;  // Left startpoint
      outf << xi << " 0 1\n";
      outf << xi << " " << ych1 << " 2\n";
      xi = -width/2.;
      outf << xi << " " << ych1 << " 2\n";
      outf << 0. << " " << ych2 << " 2\n";
      if( geom=="symxy" )
      {
          outf << "0 0 3\n";
      }
      else if( geom=="symy" )
      {
          xi = -xi;
          outf << xi << " " << ych1 << " 2\n";
          xi = length/2.;
          outf << xi << " " << ych1 << " 2\n";
          outf << xi << " 0 3\n";
      }
      //writing edges
      outf << "Edges " << nvertices
           << "\n1 2 1\n";
      for( int i=2; i<nvertices-1; ++i ){
          outf << i << " " << i+1 << " 2\n";
      }
      outf << nvertices-1 << " " << nvertices << " 3\n";
      outf << nvertices << " 1 4\n";

      outf << "hVertices\n";
      for( int i=0; i<nvertices; ++i ){
          outf << "1 ";
      }
      outf.close();

      cout << "\n\n\nCreate .bamg file " << endl;
      string command ="bamg ";
      command += " -g " + bname;
      command += " -o " + filename + ".bamg";
      char tmp[100];
      sprintf(tmp," -hmax %f ",hmax);
      command.append(tmp);
      cout << command << endl;
      exec_shell(command);

      // Create domain list file:
      string dname = filename + ".dmn";
      outf.open(dname.c_str());
      outf << "EdgeDomainNames\n"
              "4\n"
              "left\n"
              "top\n"
              "right\n"
              "bottom\n";
      outf.close();

      // Transform into geo file
      command = "bamg2geo ";
      command += filename + ".bamg ";
      command += filename + ".dmn";
      command += " > " + filename + ".geo";
      exec_shell(command);
}


void mk_mesh::bamg_squarechannel()
{
      string basename,filename,path;
      string geom;
      double length,height,amplitude,width,hmax;
      int subdivisions_top, nvertices;

      Float xi=0.,dx=0.;

      {
        path_t ppath[] = {"workdir"};
        get_from_xml(ppath, &path);
      }
      {
        path_t path[] = {"basename"};
        get_from_xml(path, &basename);
      }
      filename = path + "/" + basename;

      {
        path_t path[] = {"mesh","geom"};
        get_from_xml(path, &geom);
      }
      {
        path_t path[] = {"mesh","length"};
        get_from_xml(path, &length);
      }
      {
        path_t path[] = {"mesh","height"};
        get_from_xml(path, &height);
      }
      {
        path_t path[] = {"mesh","width"};
        get_from_xml(path,&width);
      }
      {
        path_t path[] = {"mesh","amplitude"};
        get_from_xml(path, &amplitude);
      }
      {
        path_t path[] = {"mesh","subdivisions_top"};
        get_from_xml(path, &subdivisions_top);
      }
      {
        path_t path[] = {"mesh","hmax"};
        get_from_xml(path, &hmax);
      }

      cout << "\n\t\tWorkdir              = " << path << endl
              << "\n\t\tbasename              = " << basename << endl
              << "\n\t\tmesh/geom             = " << geom << endl
              << "\n\t\tmesh/length           = " << length << endl
              << "\n\t\tmesh/height           = " << height << endl
              << "\n\t\tmesh/amplitude        = " << amplitude << endl
              << "\n\t\tmesh/hmax             = " << hmax << endl
              << "\n\t\tmesh/subdivsions_top  = " << subdivisions_top << endl;
      cout << "...done\n" << endl;

      if ("symy" == geom)
      {
              dx=length/subdivisions_top; //stepsize
              nvertices = 8;
      }
      else if ("symxy" == geom)
      {
              length = length*2;
              dx=length/(2.*subdivisions_top);
              nvertices = 6;
      }
      else
      {
              cerr << "Wrong geometry type !!!" << endl;
              std::exit(-1);
      }

      string bname = filename + "_g.msh";
      ofstream outf(bname.c_str());
      outf << "MeshVersionFormatted 0\n"
           << "AngleOfCornerBound 46\n"
           << "Dimension 2\n"
           << "Vertices " << nvertices
           << '\n';

      double ych1 = height-amplitude;
      double ych2 = height+amplitude;
      xi= -length/2.;  // Left startpoint
      outf << xi << " 0 1\n";
      outf << xi << " " << ych1 << " 2\n";
      xi = -width/2.;
      outf << xi << " " << ych1 << " 2\n";
      outf << xi << " " << ych2 << " 2\n";
      if( geom=="symxy" )
      {
          outf << "0 " << ych2 << " 2\n";
          outf << "0 0 3\n";
      }
      else if( geom=="symy" )
      {
          xi = -xi;
          outf << xi << " " << ych2 << " 2\n";
          outf << xi << " " << ych1 << " 2\n";
          xi = length/2.;
          outf << xi << " " << ych1 << " 2\n";
          outf << xi << " 0 3\n";
      }
      //writing edges
      outf << "Edges " << nvertices
           << "\n1 2 1\n";
      for( int i=2; i<nvertices-1; ++i ){
          outf << i << " " << i+1 << " 2\n";
      }
      outf << nvertices-1 << " " << nvertices << " 3\n";
      outf << nvertices << " 1 4\n";

      outf << "hVertices\n";
      for( int i=0; i<nvertices; ++i ){
          outf << "5 ";
      }
      outf.close();

      cout << "\n\n\nCreate .bamg file " << endl;
      string command ="bamg ";
      command += " -g " + bname;
      command += " -o " + filename + ".bamg";
      char tmp[100];
      sprintf(tmp," -hmax %f ",hmax);
      command.append(tmp);
      cout << command << endl;
      exec_shell(command);

      // Create domain list file:
      string dname = filename + ".dmn";
      outf.open(dname.c_str());
      outf << "EdgeDomainNames\n"
              "4\n"
              "left\n"
              "top\n"
              "right\n"
              "bottom\n";
      outf.close();

      // Transform into geo file
      command = "bamg2geo ";
      command += filename + ".bamg ";
      command += filename + ".dmn";
      command += " > " + filename + ".geo";
      exec_shell(command);
}

/**
 * \fn mk_mesh::bamg_wavychannel(XMLParser & xmlparam )
 * @brief Creates the geometry for the wavychannel from the XML file
 * @param xmlparam XML input file
 * @return flag
 */
int mk_mesh::bamg_wavychannel( )
{
    string basename,filename,path;
	string geom;
	double length,height,amplitude,hmax;
	int subdivisions_top;
	

	Float xi=0.,dx=0.;
	// ----------------------------------------------
	// Read xml file
	// ----------------------------------------------
	cout << "Parse XML file ..." << endl;

//	xmlparam.get_childelement_stream("basename") >> basename;
	//added by ali 29 Jun 2011
	{
	  path_t ppath[] = {"workdir"};
	  get_from_xml(ppath, &path);
	}
	{
	  path_t path[] = {"basename"};
	  get_from_xml(path, &basename);
	}
	filename = path + "/" + basename;	

//	xmlparam.get_childelement_stream("mesh/hmax") >> hmax;
	//added by ali 29 Jun 2011
	{
	  path_t path[] = {"mesh","geom"};
	  get_from_xml(path, &geom);
	}
        {
          path_t path[] = {"mesh","length"};
          get_from_xml(path, &length);
        }
        {
          path_t path[] = {"mesh","height"};
          get_from_xml(path, &height);
        }
        {
          path_t path[] = {"mesh","amplitude"};
          get_from_xml(path, &amplitude);
        }
        {
          path_t path[] = {"mesh","subdivisions_top"};
          get_from_xml(path, &subdivisions_top);
        }
        {
          path_t path[] = {"mesh","hmax"};
          get_from_xml(path, &hmax);
        }

	cout << "\n\t\tWorkdir              = " << path << endl
		<< "\n\t\tbasename              = " << basename << endl
		<< "\n\t\tmesh/geom             = " << geom << endl
		<< "\n\t\tmesh/length           = " << length << endl
		<< "\n\t\tmesh/height           = " << height << endl
		<< "\n\t\tmesh/amplitude        = " << amplitude << endl
		<< "\n\t\tmesh/hmax             = " << hmax << endl
		<< "\n\t\tmesh/subdivsions_top  = " << subdivisions_top << endl;
	cout << "...done" << endl<<endl;


	if ("symy" == geom)
	{
 		xi=-length/2.;	// Left startpoint
 		dx=length/subdivisions_top; //stepsize
 	}
	else if ("symxy" == geom)
	{
		length = length*2;
		xi=-length/2.;
		dx=length/(2.*subdivisions_top);
	}
	else 
	{
		cerr << "Wrong geometry type !!!" << endl;
		return -1; 
	}

	string bname = filename + "_g.msh";
	ofstream outf(bname.c_str());
	outf << "MeshVersionFormatted 0" << endl
		 << "AngleOfCornerBound 46" << endl
		 << "Dimension 2" << endl
		 << "Vertices " << itos(subdivisions_top+3)
		 << endl;
	outf << xi << " 0 " << 1 << endl; // bottom left point
	outf << xi << " " << height-amplitude << " " << 1 << endl; // top left point
	//added by ali 29 Jun 2011: change i<subdi.. to i<=sub
	//for symmy geometry use i<=subdivisions_top
	//for symxy i<subdivisons_top
	int ntop = 0;
	if( geom=="symy" )
	  ntop = subdivisions_top+1;
	else if( geom=="symxy" )
	  ntop = subdivisions_top;
	else
	  std::cerr << "problem in geometry type\n";

	for(int i=1;i<ntop;i++)
	{
	  xi=-length/2.0 + i*dx;
	  outf  << xi
			<< " "
			//<< height+0.5*amplitude*(1+cos((2*pi*xi)/length))
			<< height + 1*amplitude*cos((2*pi*xi)/length)
			<< " 2" 
			<< endl;
	}
	if (geom == "symy")
	  outf << xi << " 0 1" << endl;
 	else if (geom == "symxy" ) 
	{
		outf << "0 " << height + amplitude << " 3" << endl;
		outf << "0" << " 0 3" << endl;
	}
	
	// Write Edges:
	outf	<< "Edges "
			<< itos(subdivisions_top+3) 
			<< endl
			<< "1 2 1"
			<< endl;
	for(int i=2;i<=subdivisions_top+1;i++){
	  outf 	<< i << " "<< i+1 << " "
			<< "2" << endl;
	}
	outf << subdivisions_top+2 << " "
		 << subdivisions_top+3 <<" 3" << endl
		 << subdivisions_top+3 << " "
		 << "1 " <<" 4" << endl;
	outf << "hVertices" << endl;
	outf << "0.25 ";
	for(int i=1;i<=subdivisions_top;i++){
	  outf << length/subdivisions_top << " ";
	}
	outf << "0.25 " << "0.25 ";
	outf.close();
	
	// Create .bamg file
	cout << "\n\n\nCreate .bamg file " << endl;
	string command ="bamg ";
	command += " -g " + bname;
	command += " -o " + filename + ".bamg";
	char tmp[100];
    sprintf(tmp," -hmax %f ",hmax);
    command.append(tmp);
	cout << command << endl;
	system (command.c_str());
	
	// Create oambda file:
	cout << "\n\n\nCreate .amdba file " << endl;
	command ="bamg ";
	command += " -g " + bname;
	command += " -oamdba " + filename + ".amdba";
	command.append(tmp);
	cout << command << endl;
	system (command.c_str());
	
// 	// Create msh file:
// 	cout << "\n\n\nCreate .ambda file " << endl;
// 	command ="bamg ";
// 	command += " -g " + bname;
// 	command += " -omsh " + filename + ".msh";
// 	command.append(tmp);
// 	cout << command << endl;
// 	system (command.c_str());
	
	// Create domain list file:
	string dname = filename + ".dmn";
	outf.open(dname.c_str());
	outf << "EdgeDomainNames" << endl
		 << "4" << endl
		 << "left" << endl
		 << "top" << endl
		 << "right" << endl
		 << "bottom";
	outf.close();
	
	// Transform into geo file
	command = "bamg2geo ";
	command += filename + ".bamg ";
	command += filename + ".dmn";
	command += " > " + filename + ".geo";
	system (command.c_str());

	// Show geo
//	command = "geo -plotmtv " + filename + ".geo";
//	system(command.c_str());
	return 1; 
}

/**
 * @brief Creates the geometry for the ellipsoid from the XML file
 * @param xmlparam XML input file
 * @return flag
 */
int mk_mesh::bamg_ellipsoid( )
{
	string basename,filename,path;
	string geom;
	double area,ratio,xmax,ymax,hmax;
	int subdivisions_particle;
	

	Float thetai=0.,dtheta=0.;
	
	// ----------------------------------------------
	// Read xml file
	// ----------------------------------------------
	cout << "Parse XML file ..." << endl;

//	xmlparam.get_childelement_stream("workdir") >> path;
//	xmlparam.get_childelement_stream("basename") >> basename;
	filename = path + "/" + basename;	

//	xmlparam.get_childelement_stream("mesh/geom") >> geom;
//	xmlparam.get_childelement_stream("mesh/area") >> area;
//	xmlparam.get_childelement_stream("mesh/ratio") >> ratio;
//	xmlparam.get_childelement_stream("mesh/xmax") >> xmax;
//	xmlparam.get_childelement_stream("mesh/ymax") >> ymax;
//	xmlparam.get_childelement_stream("mesh/hmax") >> hmax;
//	xmlparam.get_childelement_stream("mesh/subdivisions_particle") >> subdivisions_particle;

	cout << "\n\t\tWorkdir              = " << path << endl
			<< "\n\t\tbasename              = " << basename << endl
			<< "\n\t\tmesh/geom             = " << geom << endl
			<< "\n\t\tmesh/area   	        = " << area << endl
			<< "\n\t\tmesh/ratio            = " << ratio << endl
			<< "\n\t\tmesh/xmax             = " << xmax << endl
			<< "\n\t\tmesh/ymax             = " << ymax << endl
			<< "\n\t\tmesh/hmax             = " << hmax << endl
			<< "\n\t\tmesh/subdivsions_particle  = " << subdivisions_particle << endl;
	cout << "...done" << endl<<endl;

	// Calculate a and b
	Float pi = 3.141592653589793238512808959406186204433;
	Float a = sqrt(area/(pi*ratio));
	Float b = area/(pi*a);

	cout 	<< "\n\t\ta              = " << a << endl
			<< "\n\t\tb              = " << b << endl
			<< "\n\t\ta*b*Pi         = " << a*b*pi << endl;
	cout << "...done" << endl<<endl;
	
	// Construct header
	string bname = filename + "_g.msh";
	ofstream outf(bname.c_str());
	outf << "MeshVersionFormatted 0" << endl
			<< "AngleOfCornerBound 46" << endl
			<< "Dimension 2" << endl
			<< "Vertices " << itos(subdivisions_particle+4)
			<< endl;
	
	if ("symxy" == geom)
	{		
		// Construct ellipsoid (label 5):
		Float dphi=pi/(2*subdivisions_particle);
		Float phi = 0;
		outf << a << " 0" << " 5" << endl;
		for(int i=1;i< subdivisions_particle;i++){
			outf  << a*cos(i*dphi)
					<< " "
					<< b*sin(i*dphi)
					<< " 5" 
					<< endl;
		}
		outf << "0 " << b << " 5" << endl;
		
		// Edge Points:
		outf << "0 " << ymax << " 1" << endl;
		outf << xmax << " " << ymax << " 2" << endl;
		outf << xmax << " 0" << " 1" << endl;
		
		// Write Edges:
		outf	<< "Edges "
				<< itos(subdivisions_particle+4) 
				<< endl;
		for(int i=1;i<=subdivisions_particle;i++)
		{
			outf 	<< i << " "<< i+1 << " "
					<< "5" << endl;
		}
		outf 	<< subdivisions_particle+1 << " " << subdivisions_particle+2 <<" 1" << endl
				<< subdivisions_particle+2 << " " << subdivisions_particle+3 << " 2" << endl
				<< subdivisions_particle+3 << " " << subdivisions_particle+4 << " 3" << endl
				<< subdivisions_particle+4 << " " << 1              << " 4" << endl;
		
		// Vertex sizes:
		outf << "hVertices" << endl;
		for(int i=0;i<=subdivisions_particle;i++){
			outf << 2*sin(dphi/2) << " ";
		}
		outf << "1 " << "1 " << "1 ";
		outf.close();
	}
	
	else 
	{
		cerr << "Wrong geometry type !!!" << endl;
		return -1; 
	}

	
	// Create .bamg file
	cout << "\n\n\nCreate .bamg file " << endl;
	string command ="bamg ";
	command += " -g " + bname;
	command += " -o " + filename + ".bamg";
	char tmp[100];
	sprintf(tmp," -hmax %f ",hmax);
	command.append(tmp);
	cout << command << endl;
	system (command.c_str());
	
	// Create oambda file:
	cout << "\n\n\nCreate .amdba file " << endl;
	command ="bamg ";
	command += " -g " + bname;
	command += " -oamdba " + filename + ".amdba";
	command.append(tmp);
	cout << command << endl;
	system (command.c_str());
	
// 	// Create msh file:
// 	cout << "\n\n\nCreate .ambda file " << endl;
// 	command ="bamg ";
// 	command += " -g " + bname;
// 	command += " -omsh " + filename + ".msh";
// 	command.append(tmp);
// 	cout << command << endl;
// 	system (command.c_str());
	
	// Create domain list file:
	string dname = filename + ".dmn";
	outf.open(dname.c_str());
	outf << "EdgeDomainNames" << endl
			<< "5" << endl
			<< "left" << endl
			<< "top" << endl
			<< "right" << endl
			<< "bottom" << endl
			<< "particle";
	outf.close();
	
	// Transform into geo file
	command = "bamg2geo ";
	command += filename + ".bamg ";
	command += filename + ".dmn";
	command += " > " + filename + ".geo";
	system (command.c_str());

	// Show geo
	//command = "geo " + filename + ".geo";
	//system(command.c_str());
	return 1; 
}
