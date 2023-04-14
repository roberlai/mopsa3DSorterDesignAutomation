#ifndef MOPSA_SIM_SETTING_H
#define MOPSA_SIM_SETTING_H

#include <string>
#include <unordered_map>

#include <mopsa/util/reader.hpp>
#include <mopsa/util/setting.hpp>
#include <mopsa/logger/logger.hpp>
#include <mopsa/headerdef.hpp>

#define SIM_SETTING_UNINIT -99999999

namespace mopsa
{
// Class : SimSettingBase
// Abstract : common simulaiton setting between 2D and 3D
// ----------------------------------------------------------------------------
class SimSettingBase : public SettingBase
{

public:

  SimSettingBase();

/*------ parameters ------*/
  double time_resolution{1.0};

  std::vector<double> dPs;

  // maximum step to simulation
  int boundary_max_timestep{100000};

  double collision_distance_ratio_threshold{0.05};

  std::string chip_name{"default_name"};
  std::string output_folder{"."};
  std::string debug_output_folder{"."};

  bool disable_wall_effect{false};
  bool disable_chip_partition{false};

  bool show_simulation_progress{false};

  int num_threads{1};

  bool dump_debug_step_file{false};
  std::vector<int> debug_step;

  int debug_range_start{SIM_SETTING_UNINIT};
  int debug_range_stop{SIM_SETTING_UNINIT};

  std::string init_position_path{""};

  bool dump_debug_collision_timestep{false};

  bool dump_detail_trajectory{true};
/*------ parameters ------*/

public:

  void dump(std::ostream &os) override;
};

// Class : SimSetting
// Abstract : 2-D simulaiton setting
// ----------------------------------------------------------------------------
class SimSetting : public SimSettingBase
{

public:
  
  SimSetting();

/*------ parameters ------*/
  std::string design_path;
  std::string mesh_nodes_path_x;
  std::string mesh_nodes_path_y;

  double alpha{1.0};
  double beta{1.45};

  point init_position{0, 0};
  velocity init_v{0, 0};

  double boundary_x_ratio{0.95};

  bool sim_gridflow_enabled{true};

  bool adv_wall_effect_enabled{true};
  double adv_wall_effect_collision_threshold{0.01};
  double adv_wall_effect_candidate_threshold{1.0};
  double adv_wall_effect_extra_shift{0.01};
  bool debug_adv_wall_effect_consistency{false};

  bool dump_debug_matlab_correlation_file{false};
/*------ parameters ------*/

public:
  void dump(std::ostream &os) override;
};

// Class : SimSetting3d
// Abstract : 3-D simulaiton setting
// ----------------------------------------------------------------------------
class SimSetting3d : public SimSettingBase
{

public:

  SimSetting3d();

/*------ parameters ------*/
  point3d init_position;
  std::string mesh_nodes_path;

  double boundary_less_z{SIM_SETTING_UNINIT};
/*------ parameters ------*/

public:

  void dump(std::ostream &os) override;

};

}

#endif
