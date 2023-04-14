#ifndef CHIP_PARTICLE_H
#define CHIP_PARTICLE_H

#include <mopsa/headerdef.hpp>
#include <mopsa/chip/design.hpp>
#include <boost/geometry.hpp>

#include <vector>

namespace mopsa
{

//  Class : Particle
//-----------------------------------------------------------------------------
class Particle
{

public:

  Particle(const point & coord, double diameter, int num_points);

  inline const point & coord() const;

  inline double diameter() const;

  inline double radius() const;

  inline double distance(const point &) const;

  void update_coord(const point & new_coord);

  bool cover(const point & node);

  bool overlap_obstacle(const Obstacle & obstacle, 
      std::vector<point> &res, bool debug);

  void dump(std::ostream &os);

private:

  point _current_coord;

  double _diameter;
  double _radius;

  int _num_points; // # of points to describe the circle
  polygon _polygon;
};

//  Class : Particle3d
//-----------------------------------------------------------------------------
class Particle3d 
{

public:
  Particle3d(const point3d & coord, double _diameter);

  inline const point3d & coord() const;

  inline double diameter() const;

  inline double radius() const;

  inline double distance(const point3d &) const;

  void update_coord(const point3d & new_coord);

  bool cover(const point3d & node);

private:
  point3d _current_coord;

  double _diameter;
  double _radius;
};

/******************************************************************************
  Inline Function
******************************************************************************/

//  Class : Particle
//-----------------------------------------------------------------------------
inline const point & 
Particle::coord() const
{
  return _current_coord;
}

inline double 
Particle::diameter() const
{
  return _diameter;
}

inline double 
Particle::radius() const
{
  return _radius;
}

inline double 
Particle::distance(const point & coord) const
{
  return boost::geometry::distance(_current_coord, coord);
}

//  Class : Particle3d
//-----------------------------------------------------------------------------
inline const point3d & 
Particle3d::coord() const
{
  return _current_coord;
}

inline double 
Particle3d::diameter() const
{
  return _diameter;
}

inline double 
Particle3d::radius() const
{
  return _radius;
}

inline double 
Particle3d::distance(const point3d & coord) const
{
  return boost::geometry::distance(_current_coord, coord);
}

}
#endif
