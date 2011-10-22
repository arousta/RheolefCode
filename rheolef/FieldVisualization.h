/**
 * @file FieldVisualization.h
 * @date Oct 18, 2011
 * @author ali
 *
 */

#ifndef FIELDVISUALIZATION_H_
#define FIELDVISUALIZATION_H_

#include "rheolef.h"
#include "rheolef/point.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <string>

const int Dim(2);

class SegmentedLine
{
public:
  typedef rheolef::point Point;
  SegmentedLine( const Point& _beg, const Point& _end, const int _n, const rheolef::geo& Th ):
    beg(_beg)
   ,end(_end)
   ,n(_n)
  {
    assert(1<n);
    const Point d( (end-beg)/(n-1) );
    Th.init_localizer( Th.boundary() );
    Point p(beg);
    for( int i=0; i<n; ++i )
    {
        //test if point is inside mesh
        rheolef::field::size_type tmp;
        if( !Th.localize(p,tmp) )
        {
          std::cerr << "out of meshhss\n";
          continue;
        }
        points.push_back(p);
        p += d;
    }
  }

  size_t size() const
  {return points.size();}

  const Point& point( size_t i ) const
  {return points[i];}

private:
  Point beg;
  Point end;
  std::vector<Point> points;
  int n;
};


template< typename FieldType >
class Field1D
{
  SegmentedLine const* line;
  std::vector<FieldType> val;
public:
  Field1D( const SegmentedLine& _line, const rheolef::field& f )
  {create(_line,f);}

  void create( const SegmentedLine& _line, const rheolef::field& f )
  {
    line = &_line;
    val.resize( line->size() );
    for( size_t i=0; i<line->size(); ++i ){
        f.evaluate( line->point(i), val[i] );
    }
  }

  const FieldType& operator()( size_t i ) const
  {return val[i];}
};


template< typename T >
void extract_field_values_on_line(
                     const rheolef::field& f
                    ,const SegmentedLine& sline
                    ,std::vector<T>& result
                    )
{
  result.resize( sline.size() );
  for( size_t i=0; i<sline.size(); ++i ){
      f.evaluate( sline.point(i), result[i] );
  }
}


template< typename DataType >
class GnuPlotWriter {};


template<>
class GnuPlotWriter<rheolef::Float>
{
public:
//  GnuPlotWriter( const char* _name ):
//    name(_name)
//  {}
  static void write( std::ostream& o, const rheolef::Float val )
  {o<<' '; o<<val;}

  static void write_header( std::ostream& o, const std::string& name )
  {o<<" "; o<<name;}
};


template<>
class GnuPlotWriter<rheolef::point>
{
public:
  static void write( std::ostream& o, const rheolef::point& val )
  {o<<' '; val.put(o,Dim);}

  static void write_header( std::ostream& o, const std::string& name )
  {
    o<<' ';
    o<<name+"x ";
    o<<name+'y';
  }
};



template<>
class GnuPlotWriter<rheolef::tensor>
{
public:
  static void write( std::ostream& o, const rheolef::tensor& val )
  {
    o<<' ';
    o<<val(x,x); o<<' ';
    o<<val(x,y); o<<' ';
    o<<val(y,x); o<<' ';
    o<<val(y,y);
  }

  static void write_header( std::ostream& o, const std::string& name )
  {
    o<<' ';
    o<<name+"xx"; o<<' ';
    o<<name+"xy"; o<<' ';
    o<<name+"yx"; o<<' ';
    o<<name+"yy";
  }

private:
  enum {x,y};
};


template< typename Record >
void write2file(
    const char* fname
   ,const SegmentedLine& line
   ,const Record& rc
   )
{
  std::ofstream o(fname);
  o <<"# x  y ";
  rc.write_header(o);
  o<<'\n';
  for( size_t i=0; i<line.size(); ++i )
  {
      line.point(i).put(o,Dim);
      o<<' ';
      rc.write(o,i);
      o<<'\n';
  }
  o.close();
}


struct TPGammaRecord
{
  typedef GnuPlotWriter<rheolef::tensor> twriter;
  typedef GnuPlotWriter<rheolef::Float>  swriter;
  Field1D<rheolef::tensor> const* T;
  Field1D<rheolef::Float>  const* P;
  Field1D<rheolef::tensor> const* Gam;

  void write_header( std::ostream& o ) const
  {
    swriter::write_header(o,"P");
    twriter::write_header(o,"Tau");
    twriter::write_header(o,"Gamma");
  }

  void write( std::ostream& o, size_t i ) const
  {
    swriter::write( o, (*P)(i) );
    twriter::write( o, (*T)(i) );
    twriter::write( o, (*Gam)(i) );
  }

};

#endif /* FIELDVISUALIZATION_H_ */
