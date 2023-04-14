#ifndef MOPSA_SIM_SIM_FLOW_GRID_H
#define MOPSA_SIM_SIM_FLOW_GRID_H

#include <mopsa/chip/chip.hpp>
#include <mopsa/util/util_funcs.h>
#include <mopsa/logger/debug.hpp>
#include <mopsa/sim/def.hpp>

#include <cmath>
#include <iostream>
#include <optional>
#include <vector>

namespace mopsa
{

class Particle;
class Particle3d;

enum class FlowBlockNodeType
{
  flow_node,
  obstacle_node,
  obstacle,
  last
};

//  Class : FlowBlock
//  Abstract : A block storage to keep the flow/obstacle nodes id
//-----------------------------------------------------------------------------
template<class P>
class FlowBlock
{

  friend class SimFlowGrid;
  friend class SimFlowGrid3d;

public:

  static constexpr int type_size = static_cast<int>(FlowBlockNodeType::last);

  inline P lower_corner() const;

  inline P upper_corner() const;

  inline std::vector<int> const * get(FlowBlockNodeType type) const;

private:

  void _make_obstacle_unique();

  P _lower_corner;
  P _upper_corner; 
  std::vector<int> _nodes;          // nodes id in Flow
  std::vector<int> _obstacle_nodes; 
  std::vector<int> _obstacles;
};

//  Class : FlowBlockGroup
//  Abstract : A group of FlowBlock
//-----------------------------------------------------------------------------
template<class P>
class FlowBlockGroup
{
public:

  FlowBlockGroup();

  inline void append(FlowBlock<P> const * block);

  inline void init(const P & grid_pos);

  inline int get(FlowBlockNodeType type, int id) const;

  inline int size(FlowBlockNodeType type) const;

  inline void clear();

  inline const P& get_grid_pos() const;

  inline auto begin() const { return _blocks.begin(); }
  inline auto end() const { return _blocks.end(); }

  void dump_nodes(
    std::ostream &os,
    FlowBlockNodeType type, 
    std::function<void(std::ostream&,int)> func_nodeid_to_str
  );

private:

  P _grid_pos;
  std::array<std::vector<int>, FlowBlock<P>::type_size> _nodes;
  std::vector<FlowBlock<P> const *> _blocks;
};

//  Class : SimFlowGrid
//  Abstract : A 2D flow grid. Divide the design to several FlowBlock
//-----------------------------------------------------------------------------
class SimFlowGrid
{
public:

  using flow_block_group_type = flow_block_group;

  SimFlowGrid(Chip *chip, double len);

  ~SimFlowGrid();
  
  bool build();

  void 
  get_adjacent_blocks(Particle const *, std::vector<flow_block const*>&) const;

  void 
  get_adjacent_blocks(Particle const *, flow_block_group_type & group) const;

private:

  inline std::pair<int, int> _node_to_grid_xy(const point & node) const;

  inline flow_block * _get_blocks(const point & node);
  inline flow_block const * _get_blocks(const point & node) const;

  inline flow_block * _get_flow_block(int x, int y);
  inline flow_block const * _get_flow_block(int x, int y) const;

  inline bool _is_valid_grid_xy(int gridx, int gridy) const;

private:

  Chip *_chip;
  double _grid_len;

  int _num_row;
  int _num_column;
  int _num_blocks;

  flow_block *_flowgrid;
};

//  Class : SimFlowGrid3d
//  Abstract : A 3D flow grid. Divide the design to several FlowBlock
//-----------------------------------------------------------------------------
class SimFlowGrid3d
{
public:

  using flow_block_group_type = flow_block_group3d;

  SimFlowGrid3d(Chip3d *chip, double len);

  ~SimFlowGrid3d();
  
  bool build();

  void 
  get_adjacent_blocks(Particle3d const *, flow_block_group_type & group) const;

  void dump_obstacle_nodes(std::ostream &os) const;
  void dump_flow_grid(std::ostream &os) const;

  inline flow_block3d const * 
  get_flow_block(
    unsigned int x, 
    unsigned int y, 
    unsigned int z
  ) const;

  inline int num_x() const { return _num_column; }
  inline int num_y() const { return _num_row; }
  inline int num_z() const { return _num_height; }
  inline double grid_len() const { return _grid_len; }

private:

  inline flow_block3d * _get_blocks(const point3d & node);
  inline flow_block3d const * _get_blocks(const point3d & node) const;

  inline flow_block3d * _get_flow_block(int x, int y, int z);
  inline flow_block3d const * _get_flow_block(int x, int y, int z) const;

  inline std::tuple<int, int, int> 
  _node_to_grid_xyz(const point3d & node) const;

  inline bool _is_valid_grid_xyz(int gridx, int gridy, int gridz) const;

private:

  Chip3d *_chip;
  double _grid_len;

  int _num_row;
  int _num_column;
  int _num_height;
  int _num_blocks;

  int _obstacle_nodes_sum{0};
  int _flow_nodes_sum{0};

  flow_block3d *_flowgrid;
};

/******************************************************************************
  Inline Function
******************************************************************************/

//  Class : FlowBlock
//-----------------------------------------------------------------------------

template<class P>
inline P
FlowBlock<P>::lower_corner() const
{
  return _lower_corner;
}

template<class P>
inline P
FlowBlock<P>::upper_corner() const
{
  return _upper_corner;
}

template<class P>
inline std::vector<int> const *
FlowBlock<P>::get(FlowBlockNodeType type) const
{
  switch (type)
  {
    case FlowBlockNodeType::flow_node:
      return &_nodes;

    case FlowBlockNodeType::obstacle_node:
      return &_obstacle_nodes;

    case FlowBlockNodeType::obstacle:
      return &_obstacles;

    case FlowBlockNodeType::last:
      ASSERT(false);
  }

  return nullptr;
}

template<class P>
inline void 
FlowBlock<P>::_make_obstacle_unique()
{
  std::sort(_obstacles.begin(), _obstacles.end());
  _obstacles.resize(
    std::unique(_obstacles.begin(), _obstacles.end()) - _obstacles.begin()
  );
}

//  Class : FlowBlockGroup
//-----------------------------------------------------------------------------

template<class P>
FlowBlockGroup<P>::FlowBlockGroup()
{
  for(size_t i=0; i<_nodes.size(); i++) {
    _nodes[i].reserve(100000);
  }

  _grid_pos.set_to_uninit();
}

template<class P>
inline void 
FlowBlockGroup<P>::append(FlowBlock<P> const * block)
{
  for(int i=0; i<flow_block::type_size; i++) {
    auto type = static_cast<FlowBlockNodeType>(i);
    auto * data = block->get(type); 
    ASSERT(data != nullptr);

    int start_id = _nodes[i].size();
    _nodes[i].resize(start_id + data->size());

    for(size_t j=0; j<data->size(); j++) {
      _nodes[i][j + start_id] = data->at(j);
    }
  }

  _blocks.push_back(block);
}

template<class P>
inline void 
FlowBlockGroup<P>::clear()
{
  _grid_pos.set_to_uninit();

  for(int i=0; i<flow_block::type_size; i++) {
    _nodes[i].clear();
  }
}

template<class P>
inline int
FlowBlockGroup<P>::size(FlowBlockNodeType type) const
{
  return _nodes[static_cast<int>(type)].size();
}

template<class P>
inline int
FlowBlockGroup<P>::get(FlowBlockNodeType type, int id) const
{
  return _nodes[static_cast<int>(type)].at(id);
}

template<class P>
void 
FlowBlockGroup<P>::dump_nodes(
  std::ostream &os,
  FlowBlockNodeType type, 
  std::function<void(std::ostream&,int)> print_node_by_id
)
{
  os << size(type) << std::endl;
  for (int i = 0; i < size(type); i++) {
    print_node_by_id(os, get(type, i));
    os << std::endl;
  }
}

template<class P>
void
FlowBlockGroup<P>::init(const P &grid_pos)
{
  clear();

  _grid_pos = grid_pos;
}

template<class P>
inline const P& 
FlowBlockGroup<P>::get_grid_pos() const
{
  return _grid_pos;
}

//  Class : FlowGrid
//-----------------------------------------------------------------------------
inline flow_block *
SimFlowGrid::_get_blocks(const point & node)
{
  auto [gridx, gridy] = _node_to_grid_xy(node);

  return _get_flow_block(gridx, gridy);
}

inline flow_block const *
SimFlowGrid::_get_blocks(const point & node) const
{
  auto [gridx, gridy] = _node_to_grid_xy(node);

  return _get_flow_block(gridx, gridy);
}

inline flow_block *
SimFlowGrid::_get_flow_block(int gridx, int gridy)
{
  if(!_is_valid_grid_xy(gridx, gridy)) return nullptr;

  return &_flowgrid[ gridy * _num_column + gridx];
}

inline flow_block const *
SimFlowGrid::_get_flow_block(int gridx, int gridy) const
{
  if(!_is_valid_grid_xy(gridx, gridy)) return nullptr;

  return &_flowgrid[ gridy * _num_column + gridx];
}

inline std::pair<int, int> 
SimFlowGrid::_node_to_grid_xy(const point & node) const
{
  auto node_shifted = node - 
    point(_chip->design().min_x(), _chip->design().min_y());

  if(float_equal(node_shifted.x(), 0))  node_shifted.x() = 0;
  if(float_equal(node_shifted.y(), 0))  node_shifted.y() = 0;

  return {
    std::floor(node_shifted.x() / _grid_len), 
    std::floor(node_shifted.y() / _grid_len)
  };
}

inline bool 
SimFlowGrid::_is_valid_grid_xy(int gridx, int gridy) const
{
  if(gridx < 0 or gridx >= _num_column) return false;

  if(gridy < 0 or gridy >= _num_row) return false;

  return true;
}

//  Class : FlowGrid3d
//-----------------------------------------------------------------------------
inline flow_block3d* 
SimFlowGrid3d::_get_blocks(const point3d & node)
{
  auto [gridx, gridy, gridz] = _node_to_grid_xyz(node);

  return _get_flow_block(gridx, gridy, gridz);
}

inline flow_block3d const * 
SimFlowGrid3d::_get_blocks(const point3d & node) const
{
  auto [gridx, gridy, gridz] = _node_to_grid_xyz(node);

  return _get_flow_block(gridx, gridy, gridz);
}


inline flow_block3d * 
SimFlowGrid3d::_get_flow_block(int gridx, int gridy, int gridz)
{
  return &_flowgrid[ 
      gridz * (_num_column * _num_row)
    + gridy * (_num_column)
    + gridx
  ];
}

inline flow_block3d const * 
SimFlowGrid3d::get_flow_block(
  unsigned int x, 
  unsigned int y, 
  unsigned int z
) 
const
{
  if(!_is_valid_grid_xyz(x, y, z)) return nullptr;

  return _get_flow_block(x, y, z);
}

inline flow_block3d const * 
SimFlowGrid3d::_get_flow_block(int gridx, int gridy, int gridz) const
{
  return &_flowgrid[ 
      gridz * (_num_column * _num_row)
    + gridy * (_num_column)
    + gridx
  ];
}

inline std::tuple<int, int, int> 
SimFlowGrid3d::_node_to_grid_xyz(const point3d & node) const
{
  auto node_shifted = node - _chip->flow().lower_corner();

  if(float_equal(node_shifted.x(), 0))  node_shifted.x() = 0;
  if(float_equal(node_shifted.y(), 0))  node_shifted.y() = 0;
  if(float_equal(node_shifted.z(), 0))  node_shifted.z() = 0;

  return { 
    std::floor(node_shifted.x() / _grid_len),
    std::floor(node_shifted.y() / _grid_len),
    std::floor(node_shifted.z() / _grid_len),
  };
}

inline bool 
SimFlowGrid3d::_is_valid_grid_xyz(int gridx, int gridy, int gridz) const
{
  if(gridx < 0 || gridx >= _num_column) return false;

  if(gridy < 0 || gridy >= _num_row) return false;

  if(gridz < 0 || gridz >= _num_height) return false;

  return true;
}

}

#endif
