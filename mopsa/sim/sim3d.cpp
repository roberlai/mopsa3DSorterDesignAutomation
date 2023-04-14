#include <mopsa/sim/sim.hpp>
#include <mopsa/sim/sim_wall_effect.hpp>

#include <omp.h>
#include <thread>

namespace mopsa
{

//  Class : Simulate3d
//-----------------------------------------------------------------------------
Simulate3d::Simulate3d(
  Chip3d *chip, 
  SimSetting3d *setting
)
:
  _chip(chip),
  _setting(setting),
  _flowgrids(nullptr)
{
  if(!std::filesystem::exists(_setting->output_folder)) {
    std::filesystem::create_directories(_setting->output_folder);
  }
}

Simulate3d::~Simulate3d()
{
  delete _flowgrids;
  delete _debug_dump_segment;
}

bool
Simulate3d::simulate()
{
  if(!_simulate_preprocess()) return false;

  LOG(INFO) << "====== Start simulating " << _setting->dPs.size() 
    << " different diameter. ======\n";

  bool ok = true;
  for(const auto & dp : _setting->dPs) {

    double grid_len = dp;
    if(_setting->disable_chip_partition) {
      grid_len = std::max({_chip->length(), _chip->height(), _chip->width()});
    }

    if(_flowgrids) _flowgrids = new(_flowgrids) SimFlowGrid3d(_chip, grid_len);
    else _flowgrids = new SimFlowGrid3d(_chip, grid_len);

    _flowgrids->build();

    _dump_flow_grid();

    std::filesystem::path output_path = _get_output_file_path(
      ("cpp_" + _setting->chip_name + "_" 
      + to_string(dp)+ ".txt")
    );

    FILE *fp = fopen(output_path.c_str(), "w");
    fprintf(fp, "diam, dump_detail_trajectory, num_init_pos\n");
    fprintf(fp, "%lf %d %ld\n",
      dp, _setting->dump_detail_trajectory, _init_position_vec.size()
    );

    int used_thread = _is_batch_mode? _used_thread_num:1;
    #pragma omp parallel for schedule(static, 1) num_threads(used_thread)
    for(size_t init_pos_id=0; 
      init_pos_id<_init_position_vec.size(); 
      init_pos_id++
    ) 
    {
      const auto &init_pos = _init_position_vec[init_pos_id];

      Particle3d particle(init_pos, dp);
      std::vector<point3d> trajectory;

      if(!_is_batch_mode && _debug_dump_segment) {
        _debug_dump_segment->clear();
        _debug_dump_segment->set_paraticle_diam(dp);
      }
      
      LOG(INFO) << "Simulate: init position id: " << init_pos_id << '\n';

      if(!_is_batch_mode && _debug_dump_segment) {
        std::filesystem::path output_path = _get_debug_file_path(
          "debug_segment_" + to_string(particle.diameter())
          + "_" + to_string(_setting->debug_range_start) + "_" 
          + to_string(_setting->debug_range_stop) + ".txt"
        );

        _debug_dump_segment->dump(output_path);
      }

      SimStatus sim_status = _simulate_low(&particle, trajectory);
      _sim_result[init_pos_id][dp] = 
        {particle.coord(), sim_status, get_current_sim_step()};

      if(!_is_batch_mode) {
        fprintf(fp, "%ld\n", trajectory.size());
        fprintf(fp, "%.9lf %.9lf %.9lf\n", 
          trajectory.front().x(), trajectory.front().y(), trajectory.front().z()
        );
        if(_setting->dump_detail_trajectory) {
          for(size_t i=1; i<trajectory.size()-1; i++) {
            const auto & point = trajectory[i];
            fprintf(fp, "%.9lf %.9lf %.9lf\n", point.x(), point.y(), point.z());
          }
        }
        fprintf(fp, "%.9lf %.9lf %.9lf\n", 
          trajectory.back().x(), trajectory.back().y(), trajectory.back().z()
        );
      }
    }

    if(_is_batch_mode) {
      for(size_t i=0; i<_init_position_vec.size(); i++) {
        auto [ok, status] = get_sim_result(i, dp);

        const auto &init_pos = _init_position_vec[i];
        const auto &final_pos = status.final_position;

        ASSERT(ok);
        fprintf(fp, "%d\n", status.total_timestep + 1);
        fprintf(fp, "%.9lf %.9lf %.9lf\n", 
          init_pos.x(), init_pos.y(), init_pos.z()
        );
        fprintf(fp, "%.9lf %.9lf %.9lf\n", 
          final_pos.x(), final_pos.y(), final_pos.z()
        );
      }
    }

    fclose(fp);
  }

  _dump_obstacle_nodes();

  return ok;
}

std::pair<bool, SimulateResult3d>
Simulate3d::get_sim_result(
  size_t init_pos_id,
  float dp
) const
{
  if(init_pos_id >= _sim_result.size()) return {false, SimulateResult3d()};

  if(auto iter = _sim_result[init_pos_id].find(dp); 
    iter != _sim_result[init_pos_id].end()) 
  {
    return {true, iter->second};
  }

  return {false, SimulateResult3d()};
}

bool
Simulate3d::_simulate_preprocess()
{
  _preprocess_multi_thread();

  _preprocess_init_position();

  _preprocess_debug();

  return true;
}

void 
Simulate3d::_preprocess_multi_thread()
{
  int max_thread = std::thread::hardware_concurrency();

  if(max_thread < _setting->num_threads) {
    LOG(WARNING) << "Cannot use " << _setting->num_threads << " threads"
      " in this machine.\n";
  }
  _used_thread_num = std::min(_setting->num_threads, max_thread);

  LOG(INFO) << "Use " << _used_thread_num << " threads\n";

  omp_set_num_threads(_used_thread_num);

  _current_sim_steps.resize(_setting->num_threads);
}

void 
Simulate3d::_preprocess_init_position()
{
  if(_setting->init_position_path == "") {
    _init_position_vec.push_back(_setting->init_position);
  }
  else {
    _read_init_position(_setting->init_position_path);
  }

  _sim_result.resize(_init_position_vec.size());

  if(_init_position_vec.size() > 1) {
    _is_batch_mode = true;
  }
}

void
Simulate3d::_preprocess_debug()
{ 
  bool enable_debug = false;

  if(_setting->dump_debug_step_file) {
    for(const auto & step : _setting->debug_step) {
      _debug_step.insert(step);
    }

    enable_debug = true;
  }

  if(_setting->debug_range_start != SIM_SETTING_UNINIT
    && _setting->debug_range_stop != SIM_SETTING_UNINIT) 
  {
    _debug_dump_segment = new SimDumpSegmentData<point3d, velocity3d>(
      _setting->debug_range_start,
      _setting->debug_range_stop
    );

    enable_debug = true;
  }

  if(_setting->dump_debug_collision_timestep) enable_debug = true;

  if(_init_position_vec.size() > 1 && enable_debug) {

    enable_debug = false;
    LOG(INFO) << "Disable debugging in batch mode\n";
    
    if(_debug_dump_segment) {
      delete _debug_dump_segment;
      _debug_dump_segment = nullptr;
    }

    _debug_step.clear();

    _setting->dump_debug_collision_timestep = false;
  }

  if(enable_debug) {
    if(_setting->debug_output_folder != ".") {
      if(std::filesystem::exists(_setting->debug_output_folder)) {
        LOG(INFO) << "Remove debug folder: " 
          << _setting->debug_output_folder << "\n";
        std::filesystem::remove_all(_setting->debug_output_folder);
      }
      LOG(INFO) << "Create debug folder" 
        << _setting->debug_output_folder << "\n";
      std::filesystem::create_directory(_setting->debug_output_folder);
    }
  }
}

int & 
Simulate3d::get_current_sim_step()
{
  return _current_sim_steps[omp_get_thread_num()];
}

const int & 
Simulate3d::get_current_sim_step() const
{
  return _current_sim_steps[omp_get_thread_num()];
}

void
Simulate3d::_read_init_position(
  const std::filesystem::path &path
)
{
  FILE *fp = fopen(path.c_str(), "r");
  LOG(INFO) << "Reading initial position from " << path << '\n';

  int num;
  (void)!fscanf(fp, "%d\n", &num);
  _init_position_vec.resize(num);
  for(int i=0; i<num; i++) {
    float x, y, z;
    (void)!fscanf(fp, "%f %f %f\n", &x, &y, &z);
    _init_position_vec[i] = {x, y, z};
  }
  fclose(fp);

  LOG(INFO) << num << " positions are loaded.\n";
}

bool 
Simulate3d::_build_flow_grids()
{
  double max_dp = 0;
  for(const auto & val : _setting->dPs) 
    max_dp = std::max(max_dp, val);

  ASSERT(_flowgrids == nullptr);

  _flowgrids = new SimFlowGrid3d(_chip, max_dp);

  return _flowgrids->build();
}

SimStatus
Simulate3d::_simulate_low(
  Particle3d *particle,
  std::vector<point3d> &trajectory
)
{
  SimStatus sim_status;

  LOG(INFO) << "Particle starts at " << particle->coord()
    << ", diameter = " << particle->diameter() << '\n';

  trajectory.clear();
  trajectory.push_back(particle->coord());

  SimFlowGrid3d::flow_block_group_type flow_group;

  bool debug(0);
  int &cur_time_step = get_current_sim_step();
  cur_time_step = 1;

  point3d previous_coord;
  int same_coord_cnt = 0;
  velocity3d pre_vel{0,0,0};
  while(true) 
  {
    debug = _is_debug_step();
    if(debug) printf("\nStep = %d\n", cur_time_step);

    velocity3d vel = _cal_particle_velocity(particle, flow_group);
    if(cur_time_step > 1 and vel==velocity3d(0,0,0)) {
      vel = pre_vel;
    } else if (vel != velocity3d(0, 0, 0)) {
      pre_vel = vel;
    }

    if(debug) printf("veolocity: %s\n", vel.to_string().c_str());
    point3d new_coord = _apply_velocity(particle, vel);
    (void)new_coord;
    if(debug) {
      printf("After Apply velocity, partical at : %s\n", 
        new_coord.to_string().c_str()
      );
    }

    point3d after_wall = _apply_wall_effect(particle, flow_group);
    (void)after_wall;
    if(debug) {
      printf("After Apply wall effect, partical at %s\n", 
        after_wall.to_string().c_str()
      );
    }

    trajectory.push_back(particle->coord());

    if(!_is_batch_mode && _is_debug_segment_step()) {
      _debug_dump_segment->add_particle_pos(particle->coord());
    }

    bool no_move = false;
    if(cur_time_step > 0) {
      double dist = boost::geometry::distance(
        particle->coord(), previous_coord
      );
      if(dist <= 0.005) {
        same_coord_cnt++;
        no_move = true;
      }
    }

    if(!no_move) {
      previous_coord = particle->coord();
      same_coord_cnt= 1;
    }

    sim_status = _get_sim_status(particle, same_coord_cnt);

    if(trajectory.size() >= 1000) {
      double dist = boost::geometry::distance(
        particle->coord(), trajectory[ trajectory.size() - 100 ]
      );

      if(dist <= 0.001)  {
        sim_status = SimStatus::Stuck;
      }
    }

    if(sim_status != SimStatus::Simulating) {
      _show_stop_reason(sim_status);
      break;
    }

    if(_setting->show_simulation_progress) {
      printf("step = %d, coord = %s\r", cur_time_step, 
        particle->coord().to_string().c_str()
      );
      std::fflush(stdout);
    }

    if(debug) {
      std::cout << "\nPlease check" << std::endl;
      int a; std::cin >> a;
    }
    cur_time_step += 1;
  }

  LOG(INFO) << "Particle ends at " << particle->coord()
    << ", total timestep = " << cur_time_step << '\n';

  return sim_status;
}

SimStatus
Simulate3d::_get_sim_status(
  Particle3d * particle, 
  int same_coord_cnt
) const
{
  if(_setting->boundary_less_z != SIM_SETTING_UNINIT &&
     particle->coord().z() <= _setting->boundary_less_z) 
  {
    return SimStatus::ArrivalExpectedBundary;
  }

  if(get_current_sim_step() > _setting->boundary_max_timestep) {
    return SimStatus::OverTimeStep;
  }

  if(same_coord_cnt >= 100) {
    return SimStatus::Stuck;
  }

  return SimStatus::Simulating;
}

void 
Simulate3d::_show_stop_reason(SimStatus status) const
{
  switch (status) {

    case SimStatus::ArrivalExpectedBundary:
      if(_setting->show_simulation_progress) printf("\n");
      LOG(INFO) << "Stop due to z < " << _setting->boundary_less_z << '\n';
      break;

    case SimStatus::OverTimeStep:
      if(_setting->show_simulation_progress) printf("\n");
      LOG(INFO) << "Stop due to simulation step > " 
        << _setting->boundary_max_timestep << '\n';
      break;

    case SimStatus::Stuck:
      LOG(INFO) << "Stop due to partical doesn't move for 10 timesteps\n";
      break;

    default:
      break;
  };
}

velocity3d
Simulate3d::_cal_particle_velocity(
  Particle3d *particle,
  flow_block_group3d &group
)
{
  bool debug = _is_debug_step();

  _flowgrids->get_adjacent_blocks(particle, group);

  int nodes_sum = group.size(FlowBlockNodeType::flow_node);

  if(debug) {

    std::filesystem::path path = _get_debug_file_path(
      (std::to_string(get_current_sim_step()) + "_" + 
      to_string(particle->diameter()) 
      + "_vel.log")
    );

    _dump_debug_file(particle, true, false, path);

    std::cout << "Calculate velocity" << std::endl;
    std::cout << "\tflow nodes = " << nodes_sum << std::endl;
  }

  //if(_is_debug_segment_step()) {
    //for(int i=0; i<nodes_sum; i++) {
      //int node_id = group.get(FlowBlockNodeType::flow_node, i);
      //const auto & node = _chip->flow().nodes()[node_id];
      //_debug_dump_segment->add_flow_node(node.coord(), node.vel());
    //}
  //}

  if(nodes_sum == 0) return {0, 0, 0}; 

  bool use_multi_thread = !_is_batch_mode;
  int used_thread = use_multi_thread? _used_thread_num:1;

  std::vector<int> cnt_t(used_thread, 0);
  std::vector<velocity3d> vel_t(used_thread, {0,0,0});

  #pragma omp parallel for num_threads(used_thread)
  for(int i=0; i<nodes_sum; i++) {
    int tid     = use_multi_thread? omp_get_thread_num():0;
    int node_id = group.get(FlowBlockNodeType::flow_node, i);

    const auto & node = _chip->flow().nodes()[node_id];
    bool covered      = particle->cover(node.coord());

    if(covered) {
      cnt_t[tid] += 1;
      vel_t[tid] += node.vel();
    }
  }

  int cnt = 0;
  velocity3d vel{0, 0, 0};
  for(int i=0; i<used_thread; i++) {
    cnt += cnt_t[i];
    vel += vel_t[i];
  }

  if(debug) {
    std::cout << "\tcovered nodes = " << cnt << std::endl;
    std::cout << "\tVelocity = " << (vel / cnt).to_string() << std::endl;
  }

  // for(int i=0; i<nodes_sum and debug; i++) {
  //   int node_id = group.get(FlowBlockNodeType::flow_node, i);
  //   const auto & node = _chip->flow().node(node_id);
  //   double dist = boost::geometry::distance(
  //     particle->coord(), node.coord()
  //   );
  //   printf("node %d: at %s, distance = %.6lf, vel = %s\n", 
  //     node_id, 
  //     (_chip->flow().node(node_id).coord().to_string().c_str()),
  //     dist, (_chip->flow().node(node_id).vel().to_string().c_str()));
  // }
  if(cnt != 0) return vel / cnt;
  else {
    std::vector<double> 
      distance_t(used_thread, std::numeric_limits<double>::max());
    std::vector<int> node_id_t(used_thread, 0);

    #pragma omp parallel for
    for(int i=0; i<nodes_sum; i++) {

      int tid     = use_multi_thread? omp_get_thread_num():0;
      int node_id = group.get(FlowBlockNodeType::flow_node, i);
      const auto & node = _chip->flow().node(node_id);

      double dist = boost::geometry::distance(
        particle->coord(), node.coord()
      );

      if(distance_t[tid] > dist) {
        distance_t[tid] = dist;
        node_id_t[tid] = node_id;
      }
    }

    for(int i=1; i<used_thread; i++) {
      if(distance_t[i] < distance_t[0]) {
        distance_t[0] = distance_t[i];
        node_id_t[0] = node_id_t[i];
      }
    }

    return _chip->flow().node(node_id_t[0]).vel();
  }
}

point3d 
Simulate3d::_apply_velocity(
  Particle3d *particle, 
  const velocity3d &vel
)
{
  point3d new_coord = particle->coord() + vel * _setting->time_resolution;

  particle->update_coord(new_coord);

  return new_coord;
}

point3d 
Simulate3d::_apply_wall_effect(
  Particle3d * particle,
  flow_block_group3d &group
)
{
  point3d new_coord = particle->coord();
  bool debug = _is_debug_step();

  if(_setting->disable_wall_effect) {
    return new_coord;
  }

  int cur_time_step = get_current_sim_step();

  bool use_multi_thread = !_is_batch_mode;
  GeneralWallEffect
    wall_effect(
      *particle, 
      *_flowgrids, 
      _chip->flow(), 
      group,
      particle->diameter() * _setting->collision_distance_ratio_threshold,
      use_multi_thread? _used_thread_num:1,
      debug
    );

  bool ok = wall_effect.run();
  if(!ok) {
    printf("\n");
    LOG(WARNING) << "Timestep: " << cur_time_step << " "
      << wall_effect.get_status_string() << '\n';
  }

  if(debug) {

    const auto & group = wall_effect.group();
    int nodes_sum = group.size(FlowBlockNodeType::obstacle_node);

    std::cout << "Calculate wall effect" << std::endl;
    std::cout << "\tobstacles nodes = " << nodes_sum << std::endl;
    std::cout << "\tdirection = " << wall_effect.adjsut_direction() 
      << std::endl;
    std::cout << "\tdist_to_adjust = " << wall_effect.distance_adjustment() 
      << std::endl;
    std::cout << "\tCollision = " 
      << (wall_effect.has_collision()? "TRUE":"FALSE") 
      << std::endl;
    std::cout << "\told position = " << particle->coord() << std::endl;
    std::cout << "\tnew position = " << wall_effect.new_coord() << std::endl;

    std::filesystem::path path = _get_debug_file_path(
      ("debug_" + std::to_string(cur_time_step) + "_"
       + to_string(particle->diameter()) + "_obstacle.log")
    );
    _dump_debug_file(particle, false, true, path);

    std::filesystem::path path2 = _get_debug_file_path(
      ("debug_" + std::to_string(cur_time_step) + "_" 
      + to_string(particle->diameter()) + "_obstacle_adjustment.log")
    );

    std::ofstream out(path2);

    point3d direction = wall_effect.adjsut_direction();
    point3d adjustment = direction * wall_effect.distance_adjustment();

    out << adjustment.x() << ' '
        << adjustment.y() << ' '
        << adjustment.z() << std::endl;

    out << wall_effect.new_coord().x() << ' '
        << wall_effect.new_coord().y() << ' '
        << wall_effect.new_coord().z() << std::endl;

    out.close();
  }

  if(!_is_batch_mode && _is_debug_segment_step()) {
    const auto & group = wall_effect.group();
    int nodes_sum = group.size(FlowBlockNodeType::obstacle_node);
    for(int i=0; i<nodes_sum; i++) {
      int node_id = group.get(FlowBlockNodeType::obstacle_node, i);
      const auto & node = _chip->flow().nodes()[node_id];
      _debug_dump_segment->add_obstacle_nodes(node.coord());
    }
  }

  // Update particle position
  if(wall_effect.has_collision()) {
    particle->update_coord(wall_effect.new_coord());

    if(_setting->dump_debug_collision_timestep) {
      std::filesystem::path output_path = _get_debug_file_path(
        _setting->chip_name + "_collision_timestep_" + 
        to_string(particle->diameter())
      );
      FILE *fp = fopen(output_path.c_str(), "a");
      fprintf(fp, "%d %lf\n", 
        cur_time_step,
        wall_effect.distance_adjustment()
      );
      fclose(fp);
    }
  }

  return wall_effect.new_coord();
}

void
Simulate3d::_is_on_the_particle_surface(
  Particle3d * particle, 
  const point3d & node,
  bool & on_surface,
  bool & in_particle
) const
{
  double distance = boost::geometry::distance(particle->coord(), node);
  double diff = (distance - particle->radius());
  double threshold = 
    particle->diameter() * _setting->collision_distance_ratio_threshold;

  on_surface = in_particle = false;
  if(std::fabs(diff) <= threshold) {
    on_surface = true;
  } 
  else if(diff <= 0) {
    in_particle = true;
  }
}

std::filesystem::path 
Simulate3d::_get_debug_file_path(const std::string &filename)
{
  std::filesystem::path debug_folder = _setting->debug_output_folder;
  return debug_folder / filename;
}

std::filesystem::path 
Simulate3d::_get_output_file_path(const std::string &filename)
{
  std::filesystem::path output_folder = _setting->output_folder;
  return output_folder / filename;
}

bool 
Simulate3d::_is_debug_step() const
{
  return _debug_step.find(get_current_sim_step()) != _debug_step.end();
}

bool 
Simulate3d::_is_debug_segment_step() const
{
  int cur_step = get_current_sim_step();
  return _debug_dump_segment != nullptr 
    && (_setting->debug_range_start <= cur_step)
    && (cur_step <= _setting->debug_range_stop)
    && !_is_batch_mode
  ;
}

void 
Simulate3d::_dump_obstacle_nodes()
{
  std::filesystem::path path = _get_output_file_path(
    _setting->chip_name + "_obstacle_nodes.log"
  );

  LOG(INFO) << "Dump obstacle nodes to " << path << '\n';

  std::ofstream os(path);
  _flowgrids->dump_obstacle_nodes(os);
  os.close();
}

void
Simulate3d::_dump_flow_grid()
{
  std::filesystem::path path = _get_output_file_path(
    _setting->chip_name + "_flow_grid_"
      + to_string(_flowgrids->grid_len())
      + ".log"
  );

  LOG(INFO) << "Dump flow grid to " << path << '\n';

  std::ofstream os(path);
  _flowgrids->dump_flow_grid(os);
  os.close();
  return;
}

void
Simulate3d::_dump_debug_file(
  Particle3d *particle,
  bool include_vel, 
  bool include_obstacle,
  const std::filesystem::path &output_path
) const
{
  std::ofstream out(output_path);

  out << std::fixed << std::setprecision(6);

  out << particle->coord().x() << ' ' 
      << particle->coord().y() << ' ' 
      << particle->coord().z() << std::endl;
  out << particle->diameter() << std::endl;

  flow_block_group3d group;
  _flowgrids->get_adjacent_blocks(particle, group);

  if(include_vel) {
    group.dump_nodes(out, FlowBlockNodeType::flow_node, 
      [this, particle](std::ostream &os, int id) {
        const auto & node = _chip->flow().node(id);
        os << node.x() << ' '
           << node.y() << ' '
           << node.z() << ' '
           << particle->cover(node.coord());
      }
    );
  }
  else out << "0\n";

  if(include_obstacle) {
    group.dump_nodes(out, FlowBlockNodeType::obstacle_node, 
      [this, particle](std::ostream &os, int id) {
        const auto & node = _chip->flow().node(id);
        bool on_surface, in_particle;

        _is_on_the_particle_surface(
          particle, node.coord(),
          on_surface, in_particle
        );

        os << node.x() << ' '
           << node.y() << ' '
           << node.z() << ' '
           << on_surface << ' '
           << in_particle;
      }
    );
  }
  else out << "0\n";

  // Dump grid line
  double INF = std::numeric_limits<double>::max();
  double nINF = std::numeric_limits<double>::min();
  double min_x = INF, min_y = INF, min_z = INF;
  double max_x = nINF, max_y = nINF, max_z = nINF;
  bool first = 1;
  for(const auto * block : group) {
    const auto & lower_corner = 
      block->lower_corner() + _chip->flow().lower_corner();
    const auto & upper_corner = 
      block->upper_corner() + _chip->flow().lower_corner();

    if(first) {
      min_x = lower_corner.x();
      min_y = lower_corner.y();
      min_z = lower_corner.z();

      max_x = upper_corner.x();
      max_y = upper_corner.y();
      max_z = upper_corner.z();
    }
    else {
      min_x = std::min(min_x, lower_corner.x());
      min_y = std::min(min_y, lower_corner.y());
      min_z = std::min(min_z, lower_corner.z());

      max_x = std::max(max_x, upper_corner.x());
      max_y = std::max(max_y, upper_corner.y());
      max_z = std::max(max_z, upper_corner.z());
    }
    first = false;
  }

  double grid_len = _flowgrids->grid_len();
  std::cout << "Grid_len: " << grid_len << std::endl;
  std::cout << "# of blocks: " << group.end() - group.begin() << std::endl;
  std::cout << min_x << " " << min_y << " " << min_z << std::endl;
  std::cout << max_x << " " << max_y << " " << max_z << std::endl;
  for(double x = min_x; x <= max_x + grid_len; x += grid_len) {
    for(double z = min_z; z <= max_z + grid_len; z += grid_len) {
      out  << std::min(x, max_x) << " " 
          << min_y << " " 
          << std::min(z, max_z) << " ";

      out  << std::min(x, max_x) << " " 
          << max_y << " " 
          << std::min(z, max_z) << "\n";
    }
  }

  for(double x = min_x; x <= max_x + grid_len; x += grid_len) {
    for(double y = min_y; y <= max_y + grid_len; y += grid_len) {
      out  << std::min(x, max_x) << " " 
          << std::min(y, max_y) << " "
          << min_z << " ";

      out << std::min(x, max_x) << " " 
          << std::min(y, max_y) << " "
          << max_z << "\n";
    }
  }

  for(double z = min_z; z <= max_z + grid_len; z += grid_len) {
    for(double y = min_y; y <= max_y + grid_len; y += grid_len) {
      out << min_x << " " 
          << std::min(y, max_y) << " "
          << std::min(z, max_z) << " ";

      out << max_x << " " 
          << std::min(y, max_y) << " "
          << std::min(z, max_z) << "\n";
    }
  }
  out.close();
}

void 
Simulate3d::_dump_trajectory(
  Particle3d *particle,
  const std::filesystem::path &output_path, 
  const std::vector<point3d> &trajectory
) const
{
  LOG(INFO) << "Dump result to " << output_path << "\n";
  FILE *fp = fopen(output_path.c_str(), "w");
  fprintf(fp, "%ld %f\n", trajectory.size(), particle->diameter());
  for(const auto & point : trajectory) {
    fprintf(fp, "%lf %lf %lf\n", point.x(), point.y(), point.z());
  }
  fclose(fp);
}

}
