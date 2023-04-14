#ifndef MOPSA_GEOMETRY_POLYGON_H
#define MOPSA_GEOMETRY_POLYGON_H

#include <boost/geometry/algorithms/detail/covered_by/interface.hpp>
#include <boost/geometry/algorithms/detail/within/interface.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <mopsa/geometry/point.hpp>
namespace mopsa
{

template<class T>
class Polygon 
{

public:

  using point_type   = T;
  using polygon_type = boost::geometry::model::polygon<point_type>;

  Polygon();
  Polygon(polygon_type && polygon);

  inline void push_back(const point_type & point);
  inline void reserve(size_t num);

  inline const auto & outer() const;
  inline auto size() const;
  inline const point_type& startpoint() const;
  inline const point_type& endpoint() const;

  inline const auto & data() const;

  inline std::vector<Polygon<point_type>> 
    intersection(const Polygon<point_type> &other) const;
  inline bool within(const point_type &p) const;

private:

  boost::geometry::model::polygon<point_type> _polygon;
};

/**************************************
  Inline Function
 **************************************/

template<class T>
Polygon<T>::Polygon()
{

}

template<class T>
Polygon<T>::Polygon(polygon_type && polygon)
{
  _polygon = std::move(polygon);
}

template<class T>
inline void 
Polygon<T>::push_back(const point_type & point)
{
  _polygon.outer().push_back(point);
}

template<class T>
inline void 
Polygon<T>::reserve(size_t num)
{
  _polygon.outer().reserve(num);
}

template<class T>
inline const auto & 
Polygon<T>::outer() const
{
  return _polygon.outer();
}

template<class T>
inline const typename Polygon<T>::point_type &
Polygon<T>::startpoint() const
{
  return _polygon.outer().front();
}

template<class T>
inline const typename Polygon<T>::point_type &
Polygon<T>::endpoint() const
{
  return _polygon.outer().back();
}

template<class T>
inline const auto & 
Polygon<T>::data() const
{
  return _polygon;
}

template<class T>
inline auto 
Polygon<T>::size() const
{
  return _polygon.outer().size();
}

template<class T>
inline std::vector<Polygon<T>> 
Polygon<T>::intersection(const Polygon<T> &other) const
{
  std::deque<polygon_type> res;

  boost::geometry::intersection(data(), other.data(), res);

  std::vector<Polygon<T>> res_polygon;
  for(auto & p : res) {
    res_polygon.emplace_back(std::move(p));
  }

  return res_polygon;
}

template<class T>
inline bool 
Polygon<T>::within(const point_type &p) const
{
  return boost::geometry::covered_by(p, data());
}

}

#endif
