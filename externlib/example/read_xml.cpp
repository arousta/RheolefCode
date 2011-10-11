/**
 * @file read_xml.cpp
 *
 * @date Jun 13, 2011
 * @author ali
 */

#include <iostream>

#include "SimpleXML.h"
#include "assert.h"


template< size_t depth >
void
xml_ex( int (&a)[depth] )
{
  for( int i=0; i<depth; ++i )
    printf("value of a(%i) is %i\n",i,a[i]);
}


void read_xml( char const* fname )
{

  TiXmlDocument doc(fname);
  if( !doc.LoadFile() ){
    printf("Can not load the file %s!\n", fname);
    exit(1);
  }

//  TiXmlElement const*const root = doc.RootElement();
  TiXmlElement const* elm = doc.RootElement()->FirstChildElement("parameters");

  path_t base_pt[] = {"basename"};
  std::string _basename;
  get_value(base_pt, elm, &_basename);
  std::cerr << "readed " << _basename << '\n';

  path_t path[] = {"mesh","bamg_options","hcoef"};
  double x;
  get_value( path, elm, &x );
  printf("readed value %f\n", x);

  //test for exist
  assert( xml_path_exist(path, elm) );

  path_t path2[] = {"StokesSolver","wrong"};
  assert( !xml_path_exist(path2, elm) );

  path_t path3[] = {"StokesSolver","tolerance"};
  assert( xml_path_exist(path3, elm) );

  //figure out size of array
  int p[] = {1,0,-1,4,8};
  xml_ex(p);
}

