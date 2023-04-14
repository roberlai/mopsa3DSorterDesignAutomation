#ifndef MOPSA_CHIP_FLOW_H
#define MOPSA_CHIP_FLOW_H

#include <vector>
#include <filesystem>

#include <mopsa/headerdef.hpp>
#include <mopsa/util/reader.hpp>

namespace mopsa
{

//  Class : FlowBase
//  Abstract : Load mesh nodes from .csv file 
//-----------------------------------------------------------------------------
class FlowBase: protected Reader
{

public:

  FlowBase() : Reader("%,", " \r\t\n") {}

protected:

  std::string read_comment() override;

  bool _is_comment_prefix() override;

};

//  Class : Flow
//  Abstract : 2D mesh nodes
//-----------------------------------------------------------------------------
class Flow : public FlowBase
{

public:

  using point_type = point;

  struct Node
  {
    Node() : _coord(0,0), _vel(0, 0) {}
    
    Node(point_type p, double vx, double vy) : 
      _coord(p), _vel(vx, vy) {}

    // Getter
    const auto &vel() const { return _vel; }
    const auto &coord() const { return _coord; }

    auto vx() const { return _vel.vx(); }
    auto vy() const { return _vel.vy(); }

    auto x()  const {return _coord.x(); }
    auto y()  const {return _coord.y(); }

    // Setter
    void vx(double val) { _vel.vx() = val; }
    void vy(double val) { _vel.vy() = val; }

    private:
      point_type _coord;
      velocity _vel;
  };

  bool load_flow(const std::filesystem::path &file_x, 
      const std::filesystem::path &file_y);

  inline int size() const;

  inline const auto & nodes() const;

  inline const Node & node(int i) const;

private:

  std::vector<Node> _nodes;
};

//  Class : Flow3d
//  Abstract : 3D mesh nodes
//-----------------------------------------------------------------------------
class Flow3d : protected FlowBase
{

public:

  using point_type = point3d;

  struct Node
  {
    Node() : _coord(0, 0, 0), _vel(0, 0, 0) {}
    
    Node(point_type p, velocity3d vel) : 
      _coord(p), _vel(vel) {}

    // Getter
    const auto &vel() const { return _vel; }
    const auto &coord() const { return _coord; }

    auto vx() const { return _vel.vx(); }
    auto vy() const { return _vel.vy(); }
    auto vz() const { return _vel.vz(); }

    auto x()  const {return _coord.x(); }
    auto y()  const {return _coord.y(); }
    auto z()  const {return _coord.z(); }

    // Setter
    void vx(double val) { _vel.vx() = val; }
    void vy(double val) { _vel.vy() = val; }
    void vz(double val) { _vel.vz() = val; }

    private:
      point_type _coord;
      velocity3d _vel;
  };

  bool load_flow(const std::filesystem::path &file);

  inline int size() const;

  inline const auto & nodes() const;

  inline const Node & node(int i) const;

  inline double length() const;

  inline double width() const;

  inline double height() const;

  inline point_type lower_corner() const;

  inline point_type upper_corner() const;

private:

  std::vector<Node> _nodes;

  double _length;
  double _width;
  double _height;

  point3d _lower_corner;
  point3d _upper_corner;
};

/**************************************
  Inline Function
 **************************************/

//  Class : Flow
//-----------------------------------------------------------------------------

inline int 
Flow::size() const
{
  return _nodes.size();
}

inline const auto &
Flow::nodes() const
{
  return _nodes;
}

inline const Flow::Node & 
Flow::node(int i) const
{
  return _nodes[i];
}

//  Class : Flow3d
//-----------------------------------------------------------------------------

inline int 
Flow3d::size() const
{
  return _nodes.size();
}

inline const auto &
Flow3d::nodes() const
{
  return _nodes;
}

inline const Flow3d::Node & 
Flow3d::node(int i) const
{
  return _nodes[i];
}

inline double
Flow3d::length() const
{
  return _length;
}

inline double
Flow3d::width() const
{
  return _width;
}

inline double
Flow3d::height() const
{
  return _height;
}

inline point3d 
Flow3d::lower_corner() const
{
  return _lower_corner;
}

inline point3d 
Flow3d::upper_corner() const
{
  return _upper_corner;
}

}
#endif
