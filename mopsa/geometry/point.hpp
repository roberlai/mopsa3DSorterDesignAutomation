#ifndef MOPSA_GEOMETRY_POINT_H
#define MOPSA_GEOMETRY_POINT_H

#include <mopsa/logger/debug.hpp>
#include <mopsa/util/util_funcs.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/register/point.hpp>

namespace mopsa
{

enum class PointUsageType
{
  POINT,
  VELOCITY
};

// Class: N-dim point.
// We allow clinets to perform arithmetic operations between different
// PointUsageType.
template<unsigned int N, 
         typename T, 
         PointUsageType U, 
         /* Client class */
         typename D
        >
class Point
{

  friend D;

  using value_type     = Point<N, T, U, D>;
  using value_type_ptr = Point<N, T, U, D> *;

public:

  inline bool operator == (const Point<N, T, U, D> &rhs) const;

  inline bool operator != (const Point<N, T, U, D> &rhs) const;
  
  inline bool operator < (const Point<N, T, U, D> &rhs) const;

  inline bool operator > (const Point<N, T, U, D> &rhs) const;

  template<PointUsageType U2, typename D2>
  inline D operator + (const Point<N, T, U2, D2> &rhs) const;

  template<PointUsageType U2, typename D2>
  inline D operator - (const Point<N, T, U2, D2> &rhs) const;

  inline D operator * (const T &factor) const;

  inline D operator / (const T &factor) const;

  template<PointUsageType U2, typename D2>
  inline D& operator += (const Point<N, T, U2, D2> &rhs);

  template<PointUsageType U2, typename D2>
  inline D& operator -= (const Point<N, T, U2, D2> &rhs);

  inline D& operator *= (const T &factor);

  inline D& operator /= (const T &factor);

  inline std::string to_string() const;

  inline std::string to_raw_string() const;

  inline double len() const;

  inline void set_to_init();

  inline void set_to_uninit();

  friend inline std::ostream & 
  operator << (std::ostream &os, const Point<N, T, U, D> &point)
  {
    os << point.to_string();
    return os;
  }

protected:

  template<unsigned int IDX>
  inline T data() const;

  template<unsigned int IDX>
  inline T & data();

private:

  Point();

  T _data[N];
};

/******************************************************************************
  Inline Function
******************************************************************************/

//  Class : Point
//-----------------------------------------------------------------------------

template<unsigned int N, typename T, PointUsageType U, typename D>
Point<N, T, U, D>::Point() 
{
  set_to_init();
}

template<unsigned int N, typename T, PointUsageType U, typename D>
template<unsigned int IDX>
inline T 
Point<N, T, U, D>::data() const
{
  if constexpr (IDX >= N) {
    ASSERT(false);
    return 0;
  }

  return _data[IDX];
}

template<unsigned int N, typename T, PointUsageType U, typename D>
template<unsigned int IDX>
T& 
Point<N, T, U, D>::data()
{
  if constexpr (IDX >= N) {
    ASSERT(false);
    return 0;
  }

  return _data[IDX];
}

template<unsigned int N, typename T, PointUsageType U, typename D>
bool 
Point<N, T, U, D>::operator == (const Point<N, T, U, D> &rhs) const
{
  if constexpr (N == 1) {
    return (_data[0] == rhs._data[0]);
  }
  else if constexpr (N == 2) {
    return     (_data[0] == rhs._data[0]) 
           and (_data[1] == rhs._data[1])
           ;
  }
  else if constexpr (N == 3) {
    return     (_data[0] == rhs._data[0]) 
           and (_data[1] == rhs._data[1]) 
           and (_data[2] == rhs._data[2]) 
           ;
  }
  else if constexpr (N >= 4) {
    bool ok = true;
    for(unsigned int i=0; i<N and ok; i++) {
      ok = (_data[i] == rhs._data[i]);
    }
    return ok;
  }
}

template<unsigned int N, typename T, PointUsageType U, typename D>
bool 
Point<N, T, U, D>::operator != (const Point<N, T, U, D> &rhs) const
{
  return !((*this) == rhs);
}

template<unsigned int N, typename T, PointUsageType U, typename D>
inline bool 
Point<N, T, U, D>::operator < (const Point<N, T, U, D> &rhs) const
{
  if constexpr (N == 1) {
    return (_data[0] < rhs._data[0]);
  }
  else if constexpr (N == 2) {
    if(_data[0] == rhs._data[0]) return _data[1] < rhs._data[1];
    else return _data[0] < rhs._data[0];
  }
  else if constexpr (N == 3) {
    if(_data[0] == rhs._data[0]) {
      if(_data[1] == rhs._data[1]) return _data[2] < rhs._data[2];
      else return _data[1] < rhs._data[1];
    }
    else return _data[0] < rhs._data[0];
  }
  else {
    for(unsigned int i=0; i<N; i++) {
      if(_data[i] != rhs._data[i]) {
        return _data[i] < rhs._data[i];
      }
    }
    return false;
  }
}

template<unsigned int N, typename T, PointUsageType U, typename D>
inline bool 
Point<N, T, U, D>::operator > (const Point<N, T, U, D> &rhs) const
{
  if constexpr (N == 1) {
    return (_data[0] > rhs._data[0]);
  }
  else if constexpr (N == 2) {
    if(_data[0] == rhs._data[0]) return _data[1] > rhs._data[1];
    else return _data[0] > rhs._data[0];
  }
  else if constexpr (N == 3) {
    if(_data[0] == rhs._data[0]) {
      if(_data[1] == rhs._data[1]) return _data[2] > rhs._data[2];
      else return _data[1] > rhs._data[1];
    }
    else return _data[0] > rhs._data[0];
  }
  else {
    for(unsigned int i=0; i<N; i++) {
      if(_data[i] != rhs._data[i]) {
        return _data[i] > rhs._data[i];
      }
    }
    return false;
  }
}

template<unsigned int N, typename T, PointUsageType U, typename D>
template<PointUsageType U2, typename D2>
D 
Point<N, T, U, D>::operator + (const Point<N, T, U2, D2> &rhs) const
{
  D point;

  if constexpr (N >= 1) {
    point._data[0] = _data[0] + ((value_type_ptr)(&rhs))->_data[0];
  }
  if constexpr (N >= 2) {
    point._data[1] = _data[1] + ((value_type_ptr)(&rhs))->_data[1];
  }
  if constexpr (N >= 3) {
    point._data[2] = _data[2] + ((value_type_ptr)(&rhs))->_data[2];
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      point._data[i] = _data[i] + ((value_type_ptr)(&rhs))->_data[i];
    }
  }

  return point;
}

template<unsigned int N, typename T, PointUsageType U, typename D>
template<PointUsageType U2, typename D2>
D 
Point<N, T, U, D>::operator - (const Point<N, T, U2, D2> &rhs) const
{
  D point;

  if constexpr (N >= 1) {
    point._data[0] = _data[0] - ((value_type_ptr)(&rhs))->_data[0];
  }
  if constexpr (N >= 2) {
    point._data[1] = _data[1] - ((value_type_ptr)(&rhs))->_data[1];
  }
  if constexpr (N >= 3) {
    point._data[2] = _data[2] - ((value_type_ptr)(&rhs))->_data[2];
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      point._data[i] = _data[i] - ((value_type_ptr)(&rhs))->_data[i];
    }
  }

  return point;
}

template<unsigned int N, typename T, PointUsageType U, typename D>
D 
Point<N, T, U, D>::operator * (const T &factor) const
{
  D point;

  if constexpr (N >= 1) {
    point._data[0] = _data[0] * factor;
  }
  if constexpr (N >= 2) {
    point._data[1] = _data[1] * factor;
  }
  if constexpr (N >= 3) {
    point._data[2] = _data[2] * factor;
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      point._data[i] = _data[i] * factor;
    }
  }

  return point;
}

template<unsigned int N, typename T, PointUsageType U, typename D>
D 
Point<N, T, U, D>::operator / (const T &factor) const
{
  D point;

  if constexpr (N >= 1) {
    point._data[0] = _data[0] / factor;
  }
  if constexpr (N >= 2) {
    point._data[1] = _data[1] / factor;
  }
  if constexpr (N >= 3) {
    point._data[2] = _data[2] / factor;
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      point._data[i] = _data[i] / factor;
    }
  }

  return point;
}

template<unsigned int N, typename T, PointUsageType U, typename D>
template<PointUsageType U2, typename D2>
D &
Point<N, T, U, D>::operator += (const Point<N, T, U2, D2> &rhs)
{

  if constexpr (N >= 1) {
    _data[0] += ((value_type_ptr)(&rhs))->_data[0];
  }
  if constexpr (N >= 2) {
    _data[1] += ((value_type_ptr)(&rhs))->_data[1];
  }
  if constexpr (N >= 3) {
    _data[2] += ((value_type_ptr)(&rhs))->_data[2];
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      _data[i] += ((value_type_ptr)(&rhs))->_data[i];
    }
  }

  return *((D*)this);
}

template<unsigned int N, typename T, PointUsageType U, typename D>
template<PointUsageType U2, typename D2>
D &
Point<N, T, U, D>::operator -= (const Point<N, T, U2, D2> &rhs)
{

  if constexpr (N >= 1) {
    _data[0] -= ((value_type_ptr)(&rhs))->_data[0];
  }
  if constexpr (N >= 2) {
    _data[1] -= ((value_type_ptr)(&rhs))->_data[1];
  }
  if constexpr (N >= 3) {
    _data[2] -= ((value_type_ptr)(&rhs))->_data[2];
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      _data[i] -= ((value_type_ptr)(&rhs))->_data[i];
    }
  }

  return *((D*)this);
}

template<unsigned int N, typename T, PointUsageType U, typename D>
D &
Point<N, T, U, D>::operator *= (const T &factor)
{
  if constexpr (N >= 1) {
    _data[0] *= factor;
  }
  if constexpr (N >= 2) {
    _data[1] *= factor;
  }
  if constexpr (N >= 3) {
    _data[2] *= factor;
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      _data[i] *= factor;
    }
  }

  return *((D*)this);
}

template<unsigned int N, typename T, PointUsageType U, typename D>
D &
Point<N, T, U, D>::operator /= (const T &factor)
{
  if constexpr (N >= 1) {
    _data[0] /= factor;
  }
  if constexpr (N >= 2) {
    _data[1] /= factor;
  }
  if constexpr (N >= 3) {
    _data[2] /= factor;
  }
  if constexpr (N >= 4) {
    for(unsigned int i=3; i<N; i++) {
      _data[i] /= factor;
    }
  }

  return *((D*)this);
}

template<unsigned int N, typename T, PointUsageType U, typename D>
double 
Point<N, T, U, D>::len() const
{
  double sum = 0;
  for(unsigned int i=0; i<N; i++) {
    sum += _data[i] * _data[i];
  }
  return std::sqrt(sum);
}

template<unsigned int N, typename T, PointUsageType U, typename D>
std::string 
Point<N, T, U, D>::to_string() const
{
  std::string val = "(";

  val += mopsa::to_string(_data[0]);

  for(unsigned i=1; i<N; i++) {
    val += ", ";
    val += mopsa::to_string(_data[i]);
  }

  val += ")";

  return val;
}

template<unsigned int N, typename T, PointUsageType U, typename D>
std::string 
Point<N, T, U, D>::to_raw_string() const
{
  std::string val;

  val += mopsa::to_string(_data[0]);

  for(unsigned i=1; i<N; i++) {
    val += " ";
    val += mopsa::to_string(_data[i]);
  }

  return val;
}

template<unsigned int N, typename T, PointUsageType U, typename D>
inline void
Point<N, T, U, D>::set_to_init()
{
  memset(_data, 0, sizeof(_data));
}

template<unsigned int N, typename T, PointUsageType U, typename D>
inline void
Point<N, T, U, D>::set_to_uninit()
{
  memset(_data, -1, sizeof(_data));
}

/******************************************************************************
  Client class
******************************************************************************/
//  Class : Point2d
//-----------------------------------------------------------------------------
template<typename T>
class Point2d : public Point<2, T, PointUsageType::POINT, Point2d<T>>
{

public:

  Point2d() {}

  Point2d(T _x, T _y) {
    x() = _x;
    y() = _y;
  }

  inline T x() const { return this->template data<0>(); }
  inline T y() const { return this->template data<1>(); }

  inline T& x() { return this->template data<0>(); }
  inline T& y() { return this->template data<1>(); }
};

//  Class : Velocity2d
//-----------------------------------------------------------------------------
template<typename T>
class Velocity2d: public Point<2, T, PointUsageType::VELOCITY, Velocity2d<T>>
{

public:

  Velocity2d() {}

  Velocity2d(T _x, T _y) {
    vx() = _x;
    vy() = _y;
  }

  inline T vx() const { return this->template data<0>(); }
  inline T vy() const { return this->template data<1>(); }

  inline T& vx() { return this->template data<0>(); }
  inline T& vy() { return this->template data<1>(); }
};

//  Class : Point3d
//-----------------------------------------------------------------------------
template<typename T>
class Point3d : public Point<3, T, PointUsageType::POINT, Point3d<T>>
{

public:

  Point3d() {}

  Point3d(T _x, T _y, T _z) { 
    x() = _x;
    y() = _y;
    z() = _z;
  }

  inline T x() const { return this->template data<0>(); }
  inline T y() const { return this->template data<1>(); }
  inline T z() const { return this->template data<2>(); }

  inline T& x() { return this->template data<0>(); }
  inline T& y() { return this->template data<1>(); }
  inline T& z() { return this->template data<2>(); }
};

//  Class : Velocity3d
//-----------------------------------------------------------------------------
template<typename T>
class Velocity3d: public Point<3, T, PointUsageType::VELOCITY, Velocity3d<T>>
{

public:

  Velocity3d() {}

  Velocity3d(T _x, T _y, T _z) {
    vx() = _x;
    vy() = _y;
    vz() = _z;
  }

  inline T vx() const { return this->template data<0>(); }
  inline T vy() const { return this->template data<1>(); }
  inline T vz() const { return this->template data<2>(); }

  inline T& vx() { return this->template data<0>(); }
  inline T& vy() { return this->template data<1>(); }
  inline T& vz() { return this->template data<2>(); }
};


}

#endif
