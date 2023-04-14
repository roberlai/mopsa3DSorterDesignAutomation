#ifndef MOPSA_CHIP_DESIGN_H
#define MOPSA_CHIP_DESIGN_H

#include <mopsa/headerdef.hpp>
#include <mopsa/geometry/polygon.hpp>
#include <mopsa/geometry/point.hpp>

#include <filesystem>

namespace mopsa
{

class Obstacle
{

public:

  Obstacle(polygon poly);

  inline const polygon & get_polygon() const;

private:

  polygon _polygon;
};

class Design
{

public:

  bool load_design(const std::filesystem::path &path);

  inline int num_obstacles() const;

  inline double width() const;

  inline double height() const;

  inline double min_x() const;

  inline double min_y() const;

  inline double max_x() const;

  inline double max_y() const;

  inline const point & upper_left() const;

  inline const point & lower_right() const;

  inline const auto & obstacles() const;

private:

  std::vector<Obstacle> _obstacles;

  // the uppder left and lower right of design
  point _upper_left, _lower_right;

  double _min_x, _max_x, _min_y, _max_y;

  double _width, _height;
};

/**************************************
  Inline Function
 **************************************/
/******************
   Obstacle
 ******************/
inline const polygon & 
Obstacle::get_polygon() const
{
  return _polygon;
}

/******************
  Desgin
 ******************/
inline int
Design::num_obstacles() const
{
  return _obstacles.size();
}

inline double 
Design::width() const
{
  return _width;
}

inline double 
Design::height() const
{
  return _height;
}

inline double 
Design::min_x() const
{
  return _min_x;
}

inline double 
Design::min_y() const
{
  return _min_y;
}

inline double 
Design::max_x() const
{
  return _max_x;
}

inline double 
Design::max_y() const
{
  return _max_y;
}

inline const point &
Design::upper_left() const
{
  return _upper_left;
}

inline const point &
Design::lower_right() const
{
  return _lower_right;
}

inline const auto & 
Design::obstacles() const
{
  return _obstacles;
}

}

#endif
