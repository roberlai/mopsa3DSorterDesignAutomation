#ifndef MOPSA_SIM_SIM_H
#define MOPSA_SIM_SIM_H

#include <mopsa/headerdef.hpp>
#include <mopsa/chip/chip.hpp>
#include <mopsa/util/reader.hpp>
#include <mopsa/sim/particle.hpp>
#include <mopsa/sim/sim_flowgrid.hpp>
#include <mopsa/sim/sim_setting.hpp>
#include <vector>

namespace mopsa
{

template<class PointType, class VelType>
class SimDumpSegmentData;

enum class SimStatus
{
  ArrivalExpectedBundary,
  OverTimeStep,
  Stuck,
  Simulating
};

template<class P = point>
struct SimulateResult
{
  P final_position;
  SimStatus status;
  int total_timestep;
};
using SimulateResult3d = SimulateResult<point3d> ;

//  Class : Simulate
//  Abstract: main class for 2D simulation
//-----------------------------------------------------------------------------
class Simulate
{

public:

  Simulate(Chip *chip, SimSetting *setting);

  ~Simulate();

  bool simulate();

private:

  bool _simulate_preprocess();

  void _preprocess_debug();

  void _preprocess_multi_thread();

  bool _build_flow_grids();

  bool _simulate_low();

  bool _stop_simulation(Particle *particle, double sim_boundary_x);

  velocity _cal_particle_velocity(std::vector<int> *covered_nodes_id);

  point _apply_velocity(Particle *particle, const velocity & vel);

  point _apply_wall_effect(Particle *particle);

  bool _apply_wall_effect_low(
    Particle *particle,
    std::pair<point, point> &two_points,
    double &portion
  );

  bool _apply_wall_effect_low_adv(
    Particle *particle,
    point & new_coord
  );

  bool _apply_wall_effect_low_general_we(
    Particle *particle,
    point & new_coord
  );

  //double _cal_adv_wall_effect_portion(
    //Particle *particle,
    //point A,
    //point B
  //);

  point _cal_wall_effect_position_based_on_A_B_porition(
    Particle *particle,
    point A,
    point B,
    double portion
  );

  bool _is_large_shift(point new_coord, point old_coord);

  bool _is_debug_step();

  bool _is_debug_segment_step() const;

  void _dump_trajectory(const std::vector<point> &trajectory, 
    const std::filesystem::path &path);

  void
  _dump_debug_file(
    Particle *particle,
    bool include_vel, 
    bool include_obstacle,
    const std::filesystem::path &output_path
  ) const;


  std::filesystem::path 
  _get_debug_file_path(const std::string &filename);

private:

  Chip *_chip;

  SimSetting *_setting;

  Particle *_particle;

  SimFlowGrid *_flowgrids;

  bool _debug_wall_effect;

  int _current_sim_step;

  std::set<int> _debug_step;

  int _used_thread_num{1};

  SimDumpSegmentData<point, velocity> *_debug_dump_segment{nullptr};
};

//  Class : Simulate3d
//  Abstract: main class for 3D simulation
//-----------------------------------------------------------------------------
class Simulate3d
{

public:

  Simulate3d(Chip3d *chip, SimSetting3d *setting);
  ~Simulate3d();

  bool simulate();

  std::pair<bool, SimulateResult3d> get_sim_result(
    size_t init_pos_id, float dp
  ) const;
 
  inline const auto & get_initial_position() const;

private:

  bool _simulate_preprocess();

  void _preprocess_debug();

  void _preprocess_multi_thread();

  void _preprocess_init_position();

  bool _build_flow_grids();

  void _dump_obstacle_nodes();

  void _dump_flow_grid();

  SimStatus _simulate_low(
    Particle3d *particle, 
    std::vector<point3d> &trajectory
  );

  velocity3d _cal_particle_velocity(
    Particle3d *particle,
    flow_block_group3d &group
  );

  point3d _apply_velocity(Particle3d * particle, const velocity3d & vel);

  point3d _apply_wall_effect(Particle3d * particle, flow_block_group3d &group);

  SimStatus _get_sim_status(Particle3d * particle, int same_coord_cnt) const;

  void _show_stop_reason(SimStatus status) const;

  void
  _is_on_the_particle_surface(
    Particle3d * particle, 
    const point3d & node,
    bool & on_surface,
    bool & in_particle
  ) const;

  bool _is_debug_step() const;
  bool _is_debug_segment_step() const;

  void _read_init_position(const std::filesystem::path &path);

  std::filesystem::path 
  _get_debug_file_path(const std::string &filename);

  std::filesystem::path 
  _get_output_file_path(const std::string &filename);

  void
  _dump_debug_file(
    Particle3d *particle,
    bool include_vel, 
    bool include_obstacle,
    const std::filesystem::path &output_path
  ) const;

  void 
  _dump_trajectory(
    Particle3d *particle,
    const std::filesystem::path &output_path, 
    const std::vector<point3d> &trajectory
  ) const;

  int & get_current_sim_step();

  const int & get_current_sim_step() const;

private:

  Chip3d *_chip;

  SimSetting3d *_setting;

  SimFlowGrid3d *_flowgrids;

  std::vector<int> _current_sim_steps; // id is thread_id

  std::set<int> _debug_step;

  bool _is_batch_mode{false};

  SimDumpSegmentData<point3d, velocity3d> *_debug_dump_segment{nullptr};

  int _used_thread_num{1};

  std::vector<std::map<float, SimulateResult3d>> _sim_result;

  std::vector<point3d> _init_position_vec;
};

//  Class : Simulate3d
//-----------------------------------------------------------------------------
inline const auto & 
Simulate3d::get_initial_position() const
{
  return _init_position_vec;
}

//  Class : SimDumpSegmentData
//-----------------------------------------------------------------------------
template<class PointType, class VelType>
class SimDumpSegmentData
{

  using point = PointType;
  using velocity = VelType;

public:

  SimDumpSegmentData(int start_timestep, int end_timestep);

  void add_obstacle_nodes(point point);
  void add_particle_pos(point point);
  void add_flow_node(point pos, velocity vel);
  void set_paraticle_diam(int val);
  void clear();

  void dump(std::filesystem::path path);

private:

  int _start_timestep;
  int _end_timestep;
  int _particle_diam;

  std::set<point> _obstacle_sets;
  std::set<std::pair<point, velocity>> _flow_node_sets;
  std::vector<point> _particle_pos;
};

//  Class : SimDumpSegmentData
//-----------------------------------------------------------------------------
template<class PointType, class VelType>
SimDumpSegmentData<PointType, VelType>::SimDumpSegmentData(
  int start_timestep,
  int end_timestep
)
:
_start_timestep(start_timestep),
_end_timestep(end_timestep)
{

}

template<class PointType, class VelType>
void 
SimDumpSegmentData<PointType, VelType>::add_obstacle_nodes(PointType point)
{
  _obstacle_sets.insert(point);
}

template<class PointType, class VelType>
void
SimDumpSegmentData<PointType, VelType>::add_particle_pos(PointType point)
{
  _particle_pos.emplace_back(point);
}

template<class PointType, class VelType>
void 
SimDumpSegmentData<PointType, VelType>::add_flow_node(
  PointType pos, 
  VelType vel
)
{
  _flow_node_sets.insert({pos, vel});
}

template<class PointType, class VelType>
void 
SimDumpSegmentData<PointType, VelType>::set_paraticle_diam(int val)
{
  _particle_diam = val;
}

template<class PointType, class VelType>
void
SimDumpSegmentData<PointType, VelType>::clear()
{
  _obstacle_sets.clear();
  _particle_pos.clear();
}

template<class PointType, class VelType>
void 
SimDumpSegmentData<PointType, VelType>::dump(std::filesystem::path path)
{
  std::ofstream os(path);

  LOG(INFO) << "Dump debug segment data, timestep: [" << _start_timestep << "," 
    << _end_timestep << "] to " << path << '\n';
  LOG(INFO) << "# of obstacle ndoes = " << _obstacle_sets.size() << '\n';
  LOG(INFO) << "# of particle postion = " << _particle_pos.size() << '\n';
  LOG(INFO) << "# of flow ndoes = " << _flow_node_sets.size() << '\n';

  os << _obstacle_sets.size() << std::endl;
  for(const auto & point : _obstacle_sets) {
    os << point.to_raw_string() << std::endl;
  }

  os << _particle_pos.size() << " " << _particle_diam << std::endl;
  for(const auto & point : _particle_pos) {
    os << point.to_raw_string() << std::endl;
  }

  os << _flow_node_sets.size() << std::endl;
  for(const auto & [pos, vel]: _flow_node_sets) {
    os << pos.to_raw_string() << " " << vel.to_raw_string() << std::endl;;
  }

  os.close();
}

}

#endif
