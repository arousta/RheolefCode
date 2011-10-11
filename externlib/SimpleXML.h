/**
 * @file SimpleXML.h
 *
 * @date Jun 14, 2011
 * @author ali
 */

#ifndef SIMPLEXML_H_
#define SIMPLEXML_H_

#include <cstdio>
#include <cstdlib>
#include <string>

#include "TinyXML++/tinyxml.h"


inline void
Tsscanf( char const* input, double*const val )
{sscanf( input, "%lf", val );}

inline void
Tsscanf( char const* input, int*const val )
{sscanf( input, "%i", val );}

inline void
Tsscanf( char const* input, size_t*const val )
{sscanf( input, "%u", val );}

inline void
Tsscanf( char const* input, std::string *const val )
{val->assign(input);}


typedef char const*const path_t;

template< int depth >
inline bool
xml_path_exist( path_t (&path)[depth], TiXmlElement const* elm )
{
  for( int i=0; i<depth; ++i ){
      elm = elm->FirstChildElement(path[i]);
      if( 0==elm )
        return false;
  }
  return true;
}

template< int depth >
void
error_if_path_is_unvalid( path_t (&path)[depth], TiXmlElement const* elm )
{
  if( !xml_path_exist(path, elm) )
  {
      printf("Error: XML path not available!\n");
      for( int i=0; i<depth; ++i )
        printf("/%s", path[i]);
      printf("\n");
      exit(2);
  }
}


template< class T, int depth >
inline void
get_value( path_t (&path)[depth], TiXmlElement const* elm, T*const x )
{
//  error_if_path_is_unvalid(path,elm);

  for( int i=0; i<depth; ++i )
      elm = elm->FirstChildElement( path[i] );

  Tsscanf( elm->GetText(), x );
}






template <typename T, std::size_t N>
inline std::size_t size_of_array( T (&)[N] ) {
   return N;
}



#endif /* SIMPLEXML_H_ */
