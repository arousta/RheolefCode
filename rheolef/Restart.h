/**
 * @file Restart.h
 * @date Oct 21, 2011
 * @author ali
 */

#ifndef RESTART_H_
#define RESTART_H_

#include "rheolef.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ctime>


class Restart
{
public:
  static bool restart_file_exists( const std::string& fname )
  {return rheolef::file_exists(fname);}

  void set_restart_saving_period( double p )
  {saving_period = p;}

  template< typename T >
  static void save_std_vector( const std::vector<T>& v, std::ofstream& o, const std::string& name )
  {
    o << '\n';
    o << name;
    o << '\n';
    o << v.size();
    o << "\n\n";
    for( size_t i=0; i<v.size(); ++i ){
        o<<v[i];
        o<<'\n';
    }
  }

  template< typename T >
  static void load_std_vector( std::vector<T>& v, std::ifstream& in )
  {
    std::string tmp;
    in >> tmp;
    size_t n;
    in >> n;
    assert( 0<=n );
    v.resize(n);
    for( size_t i=0; i<n; ++i )
        in >> v[i];
  }

  bool time_to_save_for_restart()
  {
    std::time(&now);
    double dt = std::difftime(now,last_time);
    const bool reached = ( saving_period*3600<=dt );
    if( reached )
      last_time = now;
    return reached;
  }

  void start_time_monitor()
  {std::time(&last_time);}

private:
  std::time_t last_time;
  std::time_t now;
  double saving_period;
};


#endif /* RESTART_H_ */
