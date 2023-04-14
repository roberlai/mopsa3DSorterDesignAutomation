#ifndef MOPSA_HEADERDEF_H
#define MOPSA_HEADERDEF_H

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <mopsa/geometry/polygon.hpp>
#include <mopsa/geometry/point.hpp>

namespace mopsa
{

using point    = Point2d<double>;
using velocity = Velocity2d<double>;
using polygon  = Polygon<point>;

using point3d    = Point3d<double>;
using velocity3d = Velocity3d<double>;

}

BOOST_GEOMETRY_REGISTER_POINT_2D(
  mopsa::point, double, cs::cartesian, x(), y()
);

BOOST_GEOMETRY_REGISTER_POINT_3D(
  mopsa::point3d, double, cs::cartesian, x(), y(), z()
);

#endif
