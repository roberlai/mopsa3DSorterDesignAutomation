#include <algorithm>
#include <cmath>
#include <math.h>
#include <mopsa/sim/sim_flowgrid.hpp>
#include <mopsa/sim/particle.hpp>
#include <mopsa/chip/chip.hpp>
#include <mopsa/logger/logger.hpp>
#include <mopsa/util/util_funcs.h>

namespace mopsa
{

//  Class : FlowGrid
//-----------------------------------------------------------------------------

SimFlowGrid::SimFlowGrid(
  Chip *chip,
  double len
):
  _chip(chip),
  _grid_len(len)
{

  double width  = chip->width();
  double height = chip->height();

  _num_row    = std::ceil(height / _grid_len) + 1;
  _num_column = std::ceil(width  / _grid_len) + 1;
  _num_blocks = _num_row * _num_column;

  _flowgrid = new flow_block[_num_row * _num_column];

  for(int i=0; i<_num_row; i++) {
    for(int j=0; j<_num_column; j++) {
      auto *block = _get_flow_block(j, i);
      ASSERT(block != nullptr);
      block->_lower_corner = {j*_grid_len, i*_grid_len};
      block->_upper_corner = {(j+1)*_grid_len, (i+1)*_grid_len};
    }
  }

  LOG(INFO) << "Build SimFlogGrid: grid length = " << _grid_len
    << ", Grid Size: " << _num_column << " x " << _num_row << '\n';
}

SimFlowGrid::~SimFlowGrid()
{
  delete[] _flowgrid;
}

bool
SimFlowGrid::build()
{
  int node_id = -1;
  for(const auto & node : _chip->flow().nodes()) {
    node_id += 1;
    auto * block = _get_blocks(node.coord());
    if(!block) {

      LOG(WARNING) << node.coord() << " is out of design: "
        << _chip->design().upper_left() << " x " 
        << _chip->design().lower_right() 
      << '\n';

      ASSERT( node.x() >= _chip->design().min_x()
          &&  node.x() <= _chip->design().max_x()
          &&  node.y() >= _chip->design().min_y()
          &&  node.y() <= _chip->design().max_y()
      );

      continue;
    }

    if(float_equal(node.vx(), 0) and float_equal(node.vy(), 0)) 
    {
      block->_obstacle_nodes.push_back(node_id);
    }
    else block->_nodes.push_back(node_id);
  }

  int obstacle_id = -1;
  for(const auto & obstacle: _chip->design().obstacles()) {
    obstacle_id += 1;
    for(const auto & node: obstacle.get_polygon().outer()) {
      auto * block = _get_blocks(node);
      if(!block) {
        LOG(WARNING) << node << " is out of design "
          << _chip->design().upper_left() << " x " 
          << _chip->design().lower_right()
          << '\n';

        ASSERT( node.x() < _chip->design().min_x()
             || node.x() > _chip->design().max_x()
        );

         ASSERT( node.y() < _chip->design().min_y()
             || node.y() > _chip->design().max_y()
        );

        continue;
      }
      block->_obstacles.push_back(obstacle_id);
    }
  }

  size_t maximum_nodes = 0;
  size_t maximum_obstacle_nodes = 0;

  size_t nodes_sum = 0;
  size_t obstacle_node_sum = 0;
  for(int i=0; i<_num_blocks; i++) {

    nodes_sum	      += _flowgrid[i]._nodes.size();
    obstacle_node_sum += _flowgrid[i]._obstacle_nodes.size();

    maximum_nodes = 
      std::max(maximum_nodes, _flowgrid[i]._nodes.size());

    maximum_obstacle_nodes = 
      std::max(maximum_obstacle_nodes, _flowgrid[i]._obstacle_nodes.size());

    _flowgrid[i]._make_obstacle_unique();
  }

  std::string group_name = "FlowGrid_" + to_string(_grid_len);
  Logger::add_record(group_name, "# of rows", _num_row);
  Logger::add_record(group_name, "# of columns", _num_column);
  Logger::add_record(group_name, "# of nodes", nodes_sum);
  Logger::add_record(group_name, "# of obstacle nodes", obstacle_node_sum);
  Logger::add_record(group_name, "Maximum # of nodes in block", maximum_nodes);
  Logger::add_record(group_name, "Maximum # of obstacle  nodes in block", 
      maximum_obstacle_nodes);

  return true;
}

void 
SimFlowGrid::get_adjacent_blocks(
  Particle const * particle, 
  std::vector<flow_block const*> &blocks
) const
{
  blocks.clear();
  auto [gridx, gridy] = _node_to_grid_xy(particle->coord());

  for(int i=-1; i<=1; i++) {
    for(int j=-1; j<=1; j++) {
      if(_is_valid_grid_xy(gridx + i, gridy + j))  {
        auto *block = _get_flow_block(gridx + i, gridy + j);
        blocks.push_back(block);
      }
    }
  }
}

void
SimFlowGrid::get_adjacent_blocks(
  Particle const * particle,
  flow_block_group &group
) const
{
  auto [gridx, gridy] = _node_to_grid_xy(particle->coord());
  
  if(group.get_grid_pos() == point(gridx, gridy)) {
    return;
  }
  else {
    group.init(point(gridx, gridy));
  }

  for(int i=-1; i<=1; i++) {
    for(int j=-1; j<=1; j++) {
      if(_is_valid_grid_xy(gridx + i, gridy + j))  {
        auto *block = _get_flow_block(gridx + i, gridy + j);
        group.append(block);
      }
    }
  }
}

//  Class : SimFlowGrid3d
//-----------------------------------------------------------------------------

SimFlowGrid3d::SimFlowGrid3d(
  Chip3d *chip, 
  double len
)
  :
  _chip(chip),
  _grid_len(len)
{
  double length = chip->length();
  double width  = chip->width();
  double height = chip->height();

  _num_column = std::ceil(width  / _grid_len) + 1;
  _num_row    = std::ceil(length / _grid_len) + 1;
  _num_height  = std::ceil(height / _grid_len) + 1;

  _num_blocks = _num_row * _num_column * _num_height;

  _flowgrid = new flow_block3d[_num_blocks];

  for(int i=0; i<_num_height; i++) {
    for(int j=0; j<_num_row; j++) {
      for(int k=0; k<_num_column; k++) {
        auto *block = _get_flow_block(k, j, i);
        ASSERT(block != nullptr);
        block->_lower_corner = {k*_grid_len, j*_grid_len, i*_grid_len};
        block->_upper_corner = {(k+1)*_grid_len, (j+1)*_grid_len, 
                                (i+1)*_grid_len};
      }
    }
  }

  LOG(INFO) << "Build SimFlogGrid3d: grid length = " << _grid_len
    << ", Grid Size: " << _num_column << " x " << _num_row << " x " 
    << _num_height << '\n';

}

SimFlowGrid3d::~SimFlowGrid3d()
{
  delete[] _flowgrid;
}

bool 
SimFlowGrid3d::build()
{
  int node_id = -1;
  for(const auto & node : _chip->flow().nodes()) {
    node_id += 1;
    auto * block = _get_blocks(node.coord());
    if(!block) {

      LOG(WARNING) << node.coord() << " is out of design: "
        << _chip->flow().upper_corner() << " x " 
        << _chip->flow().lower_corner() 
      << '\n';

      ASSERT( node.x() >= _chip->min_x()
           && node.x() <= _chip->max_x()
           && node.y() >= _chip->min_y()
           && node.y() <= _chip->max_y()
           && node.z() >= _chip->min_z()
           && node.z() <= _chip->max_z()
           
      );

      continue;
    }

    if(float_equal(node.vx(), 0) && float_equal(node.vy(), 0)
       && float_equal(node.vz(), 0))
    {
      block->_obstacle_nodes.push_back(node_id);
    }
    else block->_nodes.push_back(node_id);
  }

  size_t maximum_nodes = 0;
  size_t maximum_obstacle_nodes = 0;

  _flow_nodes_sum = 0;
  _obstacle_nodes_sum = 0;
  for(int i=0; i<_num_blocks; i++) {

    _flow_nodes_sum     += _flowgrid[i]._nodes.size();
    _obstacle_nodes_sum += _flowgrid[i]._obstacle_nodes.size();

    maximum_nodes = 
      std::max(maximum_nodes, _flowgrid[i]._nodes.size());

    maximum_obstacle_nodes = 
      std::max(maximum_obstacle_nodes, _flowgrid[i]._obstacle_nodes.size());

    _flowgrid[i]._make_obstacle_unique();
  }

  std::string group_name = "FlowGrid3d_" + to_string(_grid_len);
  Logger::add_record(group_name, "# of columns", _num_column);
  Logger::add_record(group_name, "# of rows", _num_row);
  Logger::add_record(group_name, "# of height", _num_height);
  Logger::add_record(group_name, "# of flow nodes", _flow_nodes_sum);
  Logger::add_record(group_name, "# of obstacle nodes", _obstacle_nodes_sum);
  Logger::add_record(group_name, "Maximum # of nodes in block", maximum_nodes);
  Logger::add_record(
    group_name,
    "Maximum # of obstacle nodes in block", 
    maximum_obstacle_nodes
  );

  return true;
}
 
void
SimFlowGrid3d::dump_obstacle_nodes(std::ostream &os) const
{
  os << _obstacle_nodes_sum << '\n';
  for(int i=0; i<_num_blocks; i++) {
    const auto * obstacle_nodes = 
      _flowgrid[i].get(FlowBlockNodeType::obstacle_node);

    for(const auto & node_id : *obstacle_nodes) {
      const auto & coord = _chip->flow().node(node_id).coord();
      os << coord.x() << " " << coord.y() << " " << coord.z() << '\n';
    }
  }
}

void
SimFlowGrid3d::dump_flow_grid(std::ostream &os) const
{
  float min_x = _chip->min_x();
  float max_x = _chip->max_x();

  float min_y = _chip->min_y();
  float max_y = _chip->max_y();

  float min_z = _chip->min_z();
  float max_z = _chip->max_z();

  for(float x = min_x; x <= max_x + _grid_len; x += _grid_len) {
    for(float z = min_z; z <= max_z + _grid_len; z += _grid_len) {
      os  << std::min(x, max_x) << " " 
          << min_y << " " 
          << std::min(z, max_z) << " ";

      os  << std::min(x, max_x) << " " 
          << max_y << " " 
          << std::min(z, max_z) << "\n";
    }
  }

  for(float x = min_x; x <= max_x + _grid_len; x += _grid_len) {
    for(float y = min_y; y <= max_y + _grid_len; y += _grid_len) {
      os  << std::min(x, max_x) << " " 
          << std::min(y, max_y) << " "
          << min_z << " ";

      os  << std::min(x, max_x) << " " 
          << std::min(y, max_y) << " "
          << max_z << "\n";
    }
  }

  for(float z = min_z; z <= max_z + _grid_len; z += _grid_len) {
    for(float y = min_y; y <= max_y + _grid_len; y += _grid_len) {
      os  << min_x << " " 
          << std::min(y, max_y) << " "
          << std::min(z, max_z) << " ";

      os  << max_x << " " 
          << std::min(y, max_y) << " "
          << std::min(z, max_z) << "\n";
    }
  }
}

void 
SimFlowGrid3d::get_adjacent_blocks(
  Particle3d const * particle, 
  flow_block_group_type & group
) const
{
  auto [gridx, gridy, gridz] = _node_to_grid_xyz(particle->coord());

  if(group.get_grid_pos() == point3d(gridx, gridy, gridz)) {
    return;
  }
  else {
    group.init( point3d(gridx, gridy, gridz) );
  }

  int len = 1;
  for(int i=-len; i<=len; i++) {
    for(int j=-len; j<=len; j++) {
      for(int k=-len; k<=len; k++) {
        if(_is_valid_grid_xyz(gridx + i, gridy + j, gridz + k))  {
          auto const *block = _get_flow_block(gridx + i, gridy + j, gridz + k);
          group.append(block);
        }
      }
    }
  }
}

}
