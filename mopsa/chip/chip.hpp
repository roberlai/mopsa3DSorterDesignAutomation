#ifndef MOPSA_CHIP_CHIP_H
#define MOPSA_CHIP_CHIP_H

#include <mopsa/headerdef.hpp>
#include <mopsa/geometry/polygon.hpp>
#include <mopsa/chip/design.hpp>
#include <mopsa/chip/flow.hpp>
#include <filesystem>

namespace mopsa
{

//  Class : Chip
//-----------------------------------------------------------------------------
class Chip
{

public:
  Chip();

  bool load_design(const std::filesystem::path &design_path);

  bool load_flow(const std::filesystem::path &flow_x_path, 
      const std::filesystem::path &flow_y_path);

  inline const Flow & flow() const;

  inline const Design & design() const;

  inline double width() const;

  inline double height() const;

private:

  Design _design;

  Flow _flow;
};

//  Class : Chip3d
//-----------------------------------------------------------------------------
class Chip3d
{

public:
  Chip3d();

  bool load_flow(const std::filesystem::path &flow_path);

  inline const Flow3d & flow() const;

  inline double width() const;

  inline double length() const;

  inline double height() const;

  inline double min_x() const;
  inline double min_y() const;
  inline double min_z() const;

  inline double max_x() const;
  inline double max_y() const;
  inline double max_z() const;

private:

  Flow3d _flow;
};


/******************************************************************************
  Inline Function
 *****************************************************************************/

//  Class : Chip
//-----------------------------------------------------------------------------
inline const Flow &
Chip::flow() const
{
  return _flow;
}

inline const Design &
Chip::design() const
{
  return _design;
}

inline double 
Chip::width() const
{
  return _design.width();
}

inline double 
Chip::height() const
{
  return _design.height();
}

//  Class : Chip3d
//-----------------------------------------------------------------------------
inline const Flow3d &
Chip3d::flow() const
{
  return _flow;
}

inline double 
Chip3d::length() const
{
  return _flow.length();
}

inline double 
Chip3d::width() const
{
  return _flow.width();
}

inline double 
Chip3d::height() const
{
  return _flow.height();
}

inline double 
Chip3d::min_x() const
{
  return _flow.lower_corner().x();
}

inline double 
Chip3d::min_y() const
{
  return _flow.lower_corner().y();
}

inline double 
Chip3d::min_z() const
{
  return _flow.lower_corner().z();
}

inline double 
Chip3d::max_x() const
{
  return _flow.upper_corner().x();
}

inline double 
Chip3d::max_y() const
{
  return _flow.upper_corner().y();
}

inline double 
Chip3d::max_z() const
{
  return _flow.upper_corner().z();
}

}

#endif

