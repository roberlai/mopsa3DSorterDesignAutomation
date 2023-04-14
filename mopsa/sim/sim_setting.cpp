#include <mopsa/sim/sim_setting.hpp>
#include <mopsa/logger/logger.hpp>

namespace mopsa
{

//  Class : SimSettingBase
//-----------------------------------------------------------------------------
SimSettingBase::SimSettingBase()
{
  _register_variable("time_resolution", 
    SettingVarType::Float, 
    time_resolution
  );

  _register_variable("dPs",
    SettingVarType::FloatVec,
    dPs
  );

  _register_variable("boundary_max_timestep", 
    SettingVarType::Int,
    boundary_max_timestep
  );

  _register_variable("collision_distance_ratio_threshold",
    SettingVarType::Float,
    collision_distance_ratio_threshold
  );

  _register_variable("chip_name",
    SettingVarType::String,
    chip_name
  );

  _register_variable("output_folder",
    SettingVarType::String,
    output_folder
  );

  _register_variable("debug_output_folder",
    SettingVarType::String,
    debug_output_folder
  );

  _register_variable("disable_wall_effect",
    SettingVarType::Boolean,
    disable_wall_effect
  );

  _register_variable("disable_chip_partition",
    SettingVarType::Boolean,
    disable_chip_partition
  );

  _register_variable("num_threads",
    SettingVarType::Int,
    num_threads
  );

  _register_variable("init_position_path",
    SettingVarType::String,
    init_position_path
  );

  _register_variable("show_simulation_progress",
    SettingVarType::Boolean,
    show_simulation_progress
  );

  _register_variable("dump_debug_step_file",
    SettingVarType::Boolean,
    dump_debug_step_file
  );

  _register_variable("debug_step",
    SettingVarType::IntVec,
    debug_step
  );

  _register_variable("debug_range_start",
    SettingVarType::Int,
    debug_range_start
  );

  _register_variable("debug_range_stop",
    SettingVarType::Int,
    debug_range_stop
  );

  _register_variable("dump_debug_collision_timestep",
    SettingVarType::Boolean,
    dump_debug_collision_timestep
  );

  _register_variable("dump_detail_trajectory",
    SettingVarType::Boolean,
   dump_detail_trajectory 
  );
}

void 
SimSettingBase::dump(std::ostream &os)
{
  os << "*************************\n";
  os << " Common Setting          \n";
  os << "*************************\n";
  os << "chip name: " << chip_name << std::endl;
  os << "output_folder: " << output_folder << std::endl;
  os << "time_resolution: " << time_resolution << std::endl;
  os << "boundary_max_timestep: " << boundary_max_timestep << std::endl;
  os << "collision_distance_ratio_threshold: " << 
    collision_distance_ratio_threshold << std::endl;
  os << "dPs: ";  _dump_vec(os, dPs); os << std::endl;
  os << "num_threads: " << num_threads << std::endl;
  os << "disable_wall_effect: " << disable_wall_effect << std::endl;
  os << "dump_debug_step_file: " << dump_debug_step_file << std::endl;
  os << "dump_debug_collision_timestep: " << 
    dump_debug_collision_timestep << std::endl;
  os << "debug_step: "; _dump_vec(os, debug_step); os << std::endl;
  os << "debug_range_start: " << debug_range_start << std::endl;
  os << "debug_range_stop: " << debug_range_stop << std::endl;
}

//  Class : SimSetting
//-----------------------------------------------------------------------------
SimSetting::SimSetting() 
{
  _register_variable("design_path",
    SettingVarType::String,
    design_path
  );

  _register_variable("mesh_nodes_path_x",
    SettingVarType::String,
    mesh_nodes_path_x
  );

  _register_variable("mesh_nodes_path_y",
    SettingVarType::String,
    mesh_nodes_path_y
  );

  _register_variable("init_position_x",
    SettingVarType::Float,
    init_position.x()
  );

  _register_variable("init_position_y",
    SettingVarType::Float,
    init_position.y()
  );

  _register_variable("start_vx",
    SettingVarType::Float,
    init_v.vx()
  );

  _register_variable("start_vy",
    SettingVarType::Float,
    init_v.vy()
  );

  _register_variable("boundary_x_ratio", 
    SettingVarType::Float,
    boundary_x_ratio
  );

  _register_variable("alpha",
    SettingVarType::Float, 
    alpha
  );
  _register_variable("beta", 
    SettingVarType::Float,
    beta
  );

  _register_variable("sim_gridflow_enabled",
    SettingVarType::Boolean,
    sim_gridflow_enabled
  );

//  _register_variable("adv_wall_effect_enabled",
//    SettingVarType::Boolean,
//    adv_wall_effect_enabled 
//  );
//  _register_variable("adv_wall_effect_collision_threshold",
//    SettingVarType::Float,
//    adv_wall_effect_collision_threshold 
//  );
//  _register_variable("adv_wall_effect_candidate_threshold",
//    SettingVarType::Float,
//    adv_wall_effect_candidate_threshold
//  );
//  _register_variable("adv_wall_effect_extra_shift",
//    SettingVarType::Float,
//    adv_wall_effect_extra_shift
//  );
//  _register_variable("debug_adv_wall_effect_consistency",
//    SettingVarType::Boolean,
//    debug_adv_wall_effect_consistency
//  );

  _register_variable("dump_debug_matlab_correlation_file",
    SettingVarType::Boolean,
    dump_debug_matlab_correlation_file
  );
}

void 
SimSetting::dump(std::ostream &os)
{
  SimSettingBase::dump(os);
  os << "*************************\n";
  os << " 2D Setting              \n";
  os << "*************************\n";

  os << "design_path: " << design_path << std::endl;
  os << "mesh_nodes_path_x: " << mesh_nodes_path_x << std::endl;
  os << "mesh_nodes_path_y: " << mesh_nodes_path_y << std::endl;

  os << "init_position: " << init_position << std::endl;
  os << "init_v: " << init_v << std::endl;
  os << "alpha: " << alpha << std::endl;
  os << "beta: " << beta << std::endl;
  os << "boundary_x_ratio: " << boundary_x_ratio << std::endl;

  os << "sim_gridflow_enabled: " << sim_gridflow_enabled << std::endl;

  os << "adv_wall_effect_enabled: " << adv_wall_effect_enabled << std::endl;
  os << "disable_wall_effect: " << disable_wall_effect << std::endl;
  //  os << "adv_wall_effect_collision_threshold: " << 
  //    adv_wall_effect_collision_threshold << std::endl;
  //  os << "adv_wall_effect_candidate_threshold: " << 
  //    adv_wall_effect_candidate_threshold << std::endl;
  //  os << "adv_wall_effect_extra_shift: " << 
  //    adv_wall_effect_extra_shift << std::endl;
  //  os << "debug_adv_wall_effect_consistency: " << 
  //    debug_adv_wall_effect_consistency << std::endl;

  os << "*************************\n\n";
}

//  Class : SimSetting3d
//-----------------------------------------------------------------------------

SimSetting3d::SimSetting3d()
{
  _register_variable("mesh_nodes_path",
    SettingVarType::String,
    mesh_nodes_path
  );

  _register_variable("init_position_x",
    SettingVarType::Float,
    init_position.x()
  );

  _register_variable("init_position_y",
    SettingVarType::Float,
    init_position.y()
  );

  _register_variable("init_position_z",
    SettingVarType::Float,
    init_position.z()
  );

  _register_variable("boundary_less_z",
    SettingVarType::Float,
    boundary_less_z
  );

}

void 
SimSetting3d::dump(std::ostream &os)
{
  SimSettingBase::dump(os);

  os << "*************************\n";
  os << " 3D Setting              \n";
  os << "*************************\n";

  os << "mesh_nodes_path: " << mesh_nodes_path << std::endl;
  os << "init_position: " << init_position << std::endl;
  os << "boundary_less_z: " << boundary_less_z << std::endl;

  os << "*************************\n\n";
}

}
