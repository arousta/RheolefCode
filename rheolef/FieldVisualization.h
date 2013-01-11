/**
 * @file FieldVisualization.h
 * @date Oct 18, 2011
 * @author ali
 */

#ifndef FIELDVISUALIZATION_H_
#define FIELDVISUALIZATION_H_

#include "rheolef.h"
#include "rheolef/point.h"

#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>
#include <string>
#include <sstream>

const int Dim(2);

inline void exec_shell( const std::string& command )
{
  std::cout << command << '\n';
  int status = system(command.c_str());
}

class GeoGauge
{
public:
  GeoGauge(): Th(0)
  {}

  GeoGauge( const rheolef::geo& _Th )
  {set_geo(_Th);}

  void set_geo( const rheolef::geo& _Th )
  {
    Th = &_Th;
    Th->init_localizer( Th->boundary() );
  }

  bool is_point_inside_geo( const rheolef::point& p ) const
  {
    rheolef::field::size_type tmp;
    return( Th->localize(p,tmp) );
  }

private:
  rheolef::geo const* Th;
};


class SegmentedLine
{
public:
  typedef rheolef::point Point;
  SegmentedLine( const Point& _beg, const Point& _end, const int _n, const rheolef::geo& Th ):
    beg(_beg)
   ,end(_end)
   ,n(_n)
   ,gauge(Th)
  {
    create_line(_beg,_end,_n);
  }

  SegmentedLine()
  {}

  void set_geo( const rheolef::geo& Th )
  {gauge.set_geo(Th);}

  void create_line( const Point& beg, const Point& end, const int n )
  {
    assert(1<n);
    points.clear();
    const Point d( (end-beg)/(n-1) );
    dl = norm(d);
    Point p(beg);
    for( int i=0; i<n; ++i )
    {
        if( gauge.is_point_inside_geo(p) )
          points.push_back(p);
        else
        {
//          std::cout << "point (";
//          p.put(std::cout,Dim);
//          std::cout << ") out of mesh\n";
        }
        p += d;
    }
    std::cout << std::endl;
  }

  const rheolef::Float spacing() const
  {return dl;}

  size_t size() const
  {return points.size();}

  const Point& point( size_t i ) const
  {return points[i];}

  void erase( std::vector<size_t>& c )
  {
    assert( c.size()<size() );

    std::vector<bool> erase_flag( size(), false );
    for( size_t i=0; i<c.size(); ++i )
        erase_flag[ c[i] ] = true;

    std::vector<Point> tmp;
    for( size_t i=0; i<size(); ++i ){
        if( !erase_flag[i] )
          tmp.push_back( points[i] );
    }
    assert( tmp.size()==size()-c.size() );
    points = tmp;
  }

  std::vector<Point>::const_iterator begin() const
  {return points.begin();}

private:
  int n;
  Point beg;
  Point end;
  double dl;
  std::vector<Point> points;
  GeoGauge gauge;
};


template< typename FieldType >
class Field1D
{
  SegmentedLine const* line;
  std::vector<FieldType> val;

public:
  Field1D()
  {}

  template< typename Evaluator >
  Field1D( const SegmentedLine& _line, const rheolef::field& f, Evaluator& E )
  {create(_line,f,E);}

  template< typename Evaluator >
  void create( const SegmentedLine& _line, const rheolef::field& f, Evaluator& E )
  {
    line = &_line;
    val.resize( line->size() );
    for( size_t i=0; i<line->size(); ++i )
        E( line->point(i), f, val[i] );
  }

  size_t size() const
  {return val.size();}

  const FieldType& at( size_t i ) const
  {return val[i];}

  const FieldType& operator()( size_t i ) const
  {return val[i];}

  typename std::vector<FieldType>::const_iterator begin() const
  {return val.begin();}
};



struct FieldEvaluator
{
  template< typename FieldType >
  void operator()( const rheolef::point& p, const rheolef::field& f, FieldType& val ) const
  {f.evaluate(p,val);}
};


template< typename T >
struct MTinyVector
{
  const T& operator()( int i ) const
  {return d[i];}

  T& operator()( int i )
  {return d[i];}

private:
  T d[Dim];
};


struct GradEvaluator
{
  GradEvaluator( rheolef::Float _dl, const rheolef::geo& Th ):
    dl(_dl)
   ,dx( rheolef::point(_dl/2.,0.) )
   ,dy( rheolef::point(0.,_dl/2.) )
   ,gauge(Th)
  {}

  template< typename FieldType >
  void operator()( const rheolef::point& p, const rheolef::field& f, MTinyVector<FieldType>& val ) const
  {
    val(0) = diffcenter<FieldType>(f,p,dx)/dl;
    val(1) = diffcenter<FieldType>(f,p,dy)/dl;
  }

  void remove_unvalid_points_of_line( SegmentedLine& line ) const
  {
    std::vector<size_t> erase_flag;
    for( size_t i=0; i<line.size(); ++i )
    {
        const rheolef::point& p( line.point(i) );
        if(
              at_least_one_of_the_points_inside(p,dx)
                                                ||
              at_least_one_of_the_points_inside(p,dy)
          )
          {
            //std::printf("point i=%3u to be removed\n",i);
            erase_flag.push_back(i);
          }
    }
    line.erase( erase_flag );
  }

private:
  bool at_least_one_of_the_points_inside( const rheolef::point& p, const rheolef::point &dr ) const
  {return !(gauge.is_point_inside_geo(p+dr) || gauge.is_point_inside_geo(p-dr));}

  bool both_points_inside( const rheolef::point& p, const rheolef::point& dr ) const
  {return (gauge.is_point_inside_geo(p+dr) && gauge.is_point_inside_geo(p-dr));}

  template< typename FieldType >
  const FieldType& diffcenter( const rheolef::field& f, const rheolef::point& r, const rheolef::point& dr ) const
  {
    FieldType fpr, fmr, fp;
    if( both_points_inside(r,dr) )
    {
      f.evaluate(r+dr,fpr);
      f.evaluate(r-dr,fmr);
      return( fpr-fmr );
    }
    else if( gauge.is_point_inside_geo(r+dr) ){
        f.evaluate(r+dr,fpr);
        f.evaluate(r,fp);
        return( fpr-fp );
    }
    else if( gauge.is_point_inside_geo(r-dr) ){
        f.evaluate(r-dr,fmr);
        f.evaluate(r,fp);
        return( fp-fmr );
    }
    else{
        printf("Both points out of mesh....\n");
        exit(-1);
    }

  }

  rheolef::Float dl;
  rheolef::point dx;
  rheolef::point dy;
  GeoGauge gauge;
};


template< typename DataType >
class GnuPlotWriter {};


template<>
class GnuPlotWriter<rheolef::Float>
{
public:
  static void write( std::ostream& o, const rheolef::Float val )
  {o<<' '; o<<val;}

  static void write_header( std::ostream& o, const std::string& name )
  {o<<' '; o<<name;}
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


//removed because intel compiler errors on that
//template<>
template< typename FieldType >
class GnuPlotWriter<MTinyVector<FieldType> >
{
public:
  typedef GnuPlotWriter<FieldType> basewriter;

  static void write( std::ostream& o, const MTinyVector<FieldType>& val )
  {
    for( int i=0; i<Dim; ++i )
        basewriter::write(o,val(i));
  }

  static void write_header( std::ostream& o, const std::string& name )
  {
    std::string str[Dim] = { "X", "Y" };
    for( int i=0; i<Dim; ++i )
      basewriter::write_header(o,name+str[i]);
  }
};


template< typename Record >
void write2file(
    const std::string& fname
   ,const SegmentedLine& line
   ,const Record& rc
   )
{
  std::ofstream o(fname.c_str());
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



struct Postprocessing_FieldsOverLine
{
  Postprocessing_FieldsOverLine()
  {}

  typedef GnuPlotWriter<rheolef::tensor> twriter;
  typedef GnuPlotWriter<rheolef::Float>  swriter;
  typedef GnuPlotWriter<rheolef::point>  vwriter;

  Field1D<rheolef::Float>   P;
  Field1D<rheolef::point>   V;
  Field1D<rheolef::tensor>  T;
  Field1D<rheolef::tensor>  Gam;


  typedef GnuPlotWriter<MTinyVector<rheolef::Float> >  sswriter;
  typedef GnuPlotWriter<MTinyVector<rheolef::point> >  vvwriter;
  typedef GnuPlotWriter<MTinyVector<rheolef::tensor> > ttwriter;

  Field1D<MTinyVector<rheolef::Float> >   GP;
  Field1D<MTinyVector<rheolef::point> >   GV;
  Field1D<MTinyVector<rheolef::tensor> >  GT;
  Field1D<MTinyVector<rheolef::tensor> >  GGam;


  void write_header( std::ostream& o ) const
  {
    swriter::write_header(o,"P");
    vwriter::write_header(o,"V");
    twriter::write_header(o,"Tau");
    twriter::write_header(o,"Gamma");

    sswriter::write_header(o,"P");
    vvwriter::write_header(o,"V");
    ttwriter::write_header(o,"Tau");
    ttwriter::write_header(o,"Gamma");
  }

  void write( std::ostream& o, size_t i ) const
  {
    swriter::write( o, P(i) );
    vwriter::write( o, V(i) );
    twriter::write( o, T(i) );
    twriter::write( o, Gam(i) );

    sswriter::write(o, GP(i) );
    vvwriter::write(o, GV(i) );
    ttwriter::write(o, GT(i) );
    ttwriter::write(o, GGam(i) );
  }

  static void gen_gnuplot_script( const std::string& basename, const std::string& script )
  {
    const std::string& gnuplot = basename+".gnu";
    std::cout << "generating plot script... " << gnuplot;
    std::ofstream o(gnuplot.c_str());
    o << script;
    o.close();
    std::cout << " done\n";
  }

  static std::string datafile( const std::string& basename )
  {return basename+".dat";}

  static void viz1( const std::string& basename, const std::string& title )
  {
    const std::string& fname = datafile(basename);

    std::string script =
  "set term png\n"
  "set output '"+basename+".png'\n"
  "set grid\n"
  "set key out\n"
  "set mxtics\n"
  "set mytics\n"
  "set title '"+title+"'\n"
  "plot '"+fname+"' using 2:20 w l lw 2 t 'Txx,x' \\\n"
  "    ,'"+fname+"' using 2:14 w l lw 2 t 'P,x'   \n";
//  "    ,'"+fname+"' using 2:7  w l lw 2 t 'Txy'   \\\n"
//  "    ,'"+fname+"' using 2:21 w l lw 2 t 'Txy,x' \\\n"
//  "    ,'"+fname+"' using 2:25 w l lw 2 t 'Txy,y' \\\n"
//  "    ,'"+fname+"' using 2:11 w l lw 2 t 'Gammaxy'\n";

    gen_gnuplot_script(basename,script);
  }

  static void viz2( const std::string& basename, const std::string& title )
  {
    const std::string& fname = datafile(basename);

    std::string script =
  "set term png\n"
  "set output '"+basename+".png'\n"
  "set size 1,1\n"
  "set origin 0,0\n"
  "set grid\n"
  "set key out\n"
  "set mxtics\n"
  "set mytics\n"
  "set multiplot\n"
  "\n"
  "set size 1,.5\n"
  "set origin 0,.5\n"
  "set title '"+title+"'\n"
  "plot '"+fname+"' using 1:6 w l lw 2 t 'Txx' \\\n"
  "    ,'"+fname+"' using 1:3 w l lw 2 t 'P'   \\\n"
  "    ,'"+fname+"' using 1:7 w l lw 2 t 'Txy' \\\n"
  "    ,'"+fname+"' using 1:11 w l lw 2 t 'Gammaxy'\n"
  "\n"
  "unset title\n"
  "set size 1,.5\n"
  "set origin 0,0\n"
  "plot '"+fname+"' using 1:20 w l lw 2 t 'Txx,x' \\\n"
  "    ,'"+fname+"' using 1:14 w l lw 2 t 'P,x'   \\\n"
  "    ,'"+fname+"' using 1:21 w l lw 2 t 'Txy,x' \\\n"
  "    ,'"+fname+"' using 1:25 w l lw 2 t 'Txy,y'\n";

    gen_gnuplot_script(basename,script);
  }

  static void viz_normalstress( const std::string& basename, const std::string& title )
  {
    const std::string& fname = datafile(basename);
    std::string script =
  "set term png\n"
  "set output '"+basename+".png'\n"
  "set grid\n"
  "set key out\n"
  "set mxtics\n"
  "set mytics\n"
  "\n"
  "set title '"+title+"'\n"
  "plot '"+fname+"' using 1:($6-$3) w l lw 2 t 'Txx-p' \\\n"
  "    ,'"+fname+"' using 1:6 w l lw 2 t 'Txx'   \\\n"
  "    ,'"+fname+"' using 1:3 w l lw 2 t 'P'\n";
    gen_gnuplot_script(basename,script);
  }
};

#endif /* FIELDVISUALIZATION_H_ */
