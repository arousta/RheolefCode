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


class Restart
{
public:
  Restart( const std::string& fname )
  {restart_file = fname;}

  bool restart_file_exists() const
  {return rheolef::file_exists(restart_file+".mfield.gz");}

  void read_field( const std::string& mark, rheolef::field& f )
  {
    in.open(restart_file,"mfield");
    in >> std::setprecision(std::numeric_limits<rheolef::Float>::digits10);
    in >> rheolef::catchmark(mark);
    in >> f;
    in.close();
  }

  void finish()
  {in.close();}

private:
  rheolef::irheostream in;
  std::string restart_file;
};


#endif /* RESTART_H_ */
