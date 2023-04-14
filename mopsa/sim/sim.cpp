#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/geometry/algorithms/detail/distance/interface.hpp>

#include <mopsa/util/util_funcs.h>
#include <mopsa/geometry/point.hpp>
#include <mopsa/sim/sim.hpp>
#include <mopsa/sim/def.hpp>
#include <mopsa/sim/sim_wall_effect.hpp>
#include <mopsa/logger/logger.hpp>

#include <omp.h>

#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <ostream>
#include <string>
#include <cstdio>
#include <thread>

namespace mopsa
{

//  Class : Simulate
//-----------------------------------------------------------------------------
Simulate::Simulate(
  Chip *chip, 
  SimSetting * setting)
  : 
    _chip(chip)
  , _setting(setting)
  , _particle(nullptr)
  , _flowgrids(nullptr)
  , _debug_wall_effect(false)
{
  if(!std::filesystem::exists(_setting->output_folder)) {
    std::filesystem::create_directories(_setting->output_folder);
  }
}

Simulate::~Simulate()
{
  delete _flowgrids;
  delete _debug_dump_segment;
}

bool 
Simulate::_build_flow_grids()
{
  double max_dp = 0;
  for(const auto & val : _setting->dPs) max_dp = std::max(max_dp, val);

  _flowgrids = new SimFlowGrid(_chip, max_dp);

  return _flowgrids->build();
}

bool
Simulate::_simulate_preprocess()
{
  _preprocess_multi_thread();

  _preprocess_debug();
  
  return true;
}

void 
Simulate::_preprocess_multi_thread()
{
  int max_thread = std::thread::hardware_concurrency();

  if(max_thread < _setting->num_threads) {
    LOG(WARNING) << "Cannot use " << _setting->num_threads << " threads"
      " in this machine.\n";
  }
  _used_thread_num = std::min(_setting->num_threads, max_thread);

  LOG(INFO) << "Use " << _used_thread_num << " threads\n";

  omp_set_num_threads(_used_thread_num);
}

void
Simulate::_preprocess_debug()
{ 
  bool enable_debug = false;

  if(_setting->dump_debug_step_file) {
    for(const auto & step : _setting->debug_step) {
      _debug_step.insert(step);
    }

    enable_debug = true;
  }

  if(_setting->dump_debug_collision_timestep) enable_debug = true;

  if(_setting->dump_debug_matlab_correlation_file) enable_debug = true;

  if(_setting->debug_range_start != SIM_SETTING_UNINIT
    && _setting->debug_range_stop != SIM_SETTING_UNINIT) 
  {
    _debug_dump_segment = new SimDumpSegmentData<point, velocity>(
      _setting->debug_range_start,
      _setting->debug_range_stop
    );

    enable_debug = true;
  }

  if(enable_debug) {
    if(_setting->debug_output_folder != ".") {
      //if(std::filesystem::exists(_setting->debug_output_folder)) {
        //LOG(INFO) << "Remove " << _setting->debug_output_folder << "\n";
        //std::filesystem::remove_all(_setting->debug_output_folder);
      //}
      //LOG(INFO) << "Create " << _setting->debug_output_folder << "\n";
      std::filesystem::create_directory(_setting->debug_output_folder);
    }
  }
}

bool 
Simulate::simulate()
{
  LOG(INFO) << "====== Start simulating " << _setting->dPs.size() 
    << " different diameter. ======\n";

  if(!_simulate_preprocess()) return false;

  if(_chip->flow().size() == 0) {
    LOG(ERROR) << "No flow is loaded\n";
    return false;
  }

  bool ok = true;
  for(const auto &dp : _setting->dPs) {

    if(_setting->sim_gridflow_enabled) {
      if(_flowgrids) _flowgrids = new(_flowgrids) SimFlowGrid(_chip, dp);
      else _flowgrids = new SimFlowGrid(_chip, dp);

      _flowgrids->build();
    }

    Particle particle(_setting->init_position, dp, 200);

    if(_debug_dump_segment) {
      _debug_dump_segment->clear();
      _debug_dump_segment->set_paraticle_diam(dp);
    }

    _particle = &particle;
    if(!_simulate_low()) {
      ok = false;
      break;
    }

    if(_debug_dump_segment) {
      std::filesystem::path output_path = _get_debug_file_path(
        "debug_segment_" + to_string(_particle->diameter())
        + "_" + to_string(_setting->debug_range_start) + "_" 
        + to_string(_setting->debug_range_stop) + ".txt"
      );

      _debug_dump_segment->dump(output_path);
    }

    _particle = nullptr;
  }

  return ok;
}

bool 
Simulate::_simulate_low()
{
  std::vector<point> trajectory;
  double sim_boundary_x = 0;
  velocity pre_vel;
  bool dump_matlab_correlation = _setting->dump_debug_matlab_correlation_file;
  FILE* debug_fp = nullptr;
  std::vector<int> *covered_nodes_id = nullptr;

  LOG(INFO) << "Particle starts at " << _particle->coord()
    << ", diameter = " << _particle->diameter() << '\n';

  trajectory.push_back(_particle->coord());

  sim_boundary_x = _setting->boundary_x_ratio * _chip->design().width();

  pre_vel = _setting->init_v;

  if(dump_matlab_correlation) {
    char name[100];
    sprintf(name, "debug_%s_%s.log", _setting->chip_name.c_str(), 
      to_string(_particle->diameter()).c_str());
    const auto & path = _get_debug_file_path(name);
    debug_fp = fopen(path.c_str(), "w");
    covered_nodes_id = new std::vector<int>;
    LOG(INFO) << "Dump matlab correlation file to " << path << '\n';
  }

  _current_sim_step = 0;
  while(true) 
  {
    _current_sim_step += 1;
    bool debug = _is_debug_step();

    point cur_pos = _particle->coord();

    if(debug) printf("\n");
    _debug_wall_effect = debug;

    velocity vel = _cal_particle_velocity(covered_nodes_id);

    vel.vx() = vel.vx() * _setting->alpha;
    vel.vy() = vel.vy() * _setting->beta;

    if(float_equal(vel.vx(), 0) and float_equal(vel.vy(), 0)) {
      vel = pre_vel;
    }

    pre_vel = vel;

    point next_pos = _apply_velocity(_particle, vel);
    (void)next_pos;

    point after_wall = _apply_wall_effect(_particle); 
    (void)after_wall;

    if(dump_matlab_correlation) {
      //std::sort(covered_nodes_id->begin(), covered_nodes_id->end());
      fprintf(debug_fp, "TimeStep: %d\n", _current_sim_step);
      fprintf(debug_fp, "Current position: %f %f\n", cur_pos.x(), cur_pos.y());
      fprintf(debug_fp, "Current velocity: %f %f\n", vel.vx(), vel.vy());
      fprintf(debug_fp, "# of covered nodes: %lu\n", covered_nodes_id->size());
      fprintf(debug_fp, "Next position: %f %f\n", next_pos.x(), next_pos.y());
      fprintf(debug_fp, "After wall effect: %f %f\n", 
        after_wall.x(), after_wall.y()
      );
    }

    trajectory.push_back(_particle->coord());
    if(_is_debug_segment_step()) {
      _debug_dump_segment->add_particle_pos(_particle->coord());
    }
    
    if(_stop_simulation(_particle, sim_boundary_x)) break;

    if(_setting->show_simulation_progress) {
      printf("step = %d, coord = %s\r", _current_sim_step, 
          _particle->coord().to_string().c_str());
      std::fflush(stdout);
    }

    if(debug) { 
      std::cout << "\nPlease check" << std::endl;
      int a; std::cin >> a;
    }
  }

  if(dump_matlab_correlation) {
    fclose(debug_fp);
    delete covered_nodes_id;
  }


  LOG(INFO) << "Particle ends at " << _particle->coord() << " "
    << " diameter = " << _particle->diameter() << '\n';

  // Dump output
  std::filesystem::path output_path = 
      std::filesystem::path(_setting->output_folder) 
    / ("cpp_" + _setting->chip_name + "_" 
       + to_string(_particle->diameter())+ ".txt");

  _dump_trajectory(trajectory, output_path);
  printf("\n");

  return true;
}

bool 
Simulate::_stop_simulation(Particle *particle, double sim_boundary_x)
{
  if(_current_sim_step > _setting->boundary_max_timestep) {
    if(_setting->show_simulation_progress) printf("\n");
    LOG(INFO) << "Stop due to simulation step > " 
    << _setting->boundary_max_timestep << '\n';
    return true;
  }

  if(_particle->coord().x() >= sim_boundary_x) {
    if(_setting->show_simulation_progress) printf("\n");
    LOG(INFO) << "Stop due to x >= bondary: " 
      << sim_boundary_x << '\n';
    return true;
  }

  return false;
}

velocity 
Simulate::_cal_particle_velocity(std::vector<int> *covered_nodes_id)
{
  velocity vel(0, 0);
  int cnt = 0;
  flow_block_group block_group;
  int nodes_sum = 0;
  bool debug = _is_debug_step();

  if(covered_nodes_id) covered_nodes_id->clear();

  if(_setting->sim_gridflow_enabled) {

    _flowgrids->get_adjacent_blocks(_particle, block_group);

    nodes_sum = block_group.size(FlowBlockNodeType::flow_node);

    std::vector<int> cnt_t(_used_thread_num, 0);
    std::vector<velocity> vel_t(_used_thread_num, {0,0});

    if(debug) {
      printf("%s: nodes_sum = %d\n", __FUNCTION__, nodes_sum);
    }
    int used_thread_num = _used_thread_num;
    // TODO
    if(nodes_sum <= 400) used_thread_num = 1;
    //else if(nodes_sum <= 1000) used_thread_num = 2;
    //else if(nodes_sum <= 2000) used_thread_num = 4;
    //else if(nodes_sum <= 4000) used_thread_num = 8;
    //else if(nodes_sum <= 8000) used_thread_num = 16;

    #pragma omp parallel for num_threads(used_thread_num)
    for(int i=0; i<nodes_sum; i++) {
      int tid      = omp_get_thread_num();
      int node_id  = block_group.get(FlowBlockNodeType::flow_node, i);

      const auto & node = _chip->flow().nodes()[node_id];
      bool covered      = _particle->cover(node.coord());

      if(covered) {
        if(covered_nodes_id) covered_nodes_id->push_back(i);
        cnt_t[tid] += 1;
        vel_t[tid] += node.vel();
      }
    }

    for(int i=0; i<used_thread_num; i++) {
      cnt += cnt_t[i];
      vel += vel_t[i];
    }
  }
  else {
    #pragma omp parallel for
    for(size_t i = 0; i<_chip->flow().nodes().size(); i++) {
      const auto & node = _chip->flow().nodes()[i];
      bool covered = _particle->cover(node.coord());

      if(covered) {
        #pragma omp critical
        {
          if(covered_nodes_id) covered_nodes_id->push_back(i);
          cnt += 1;
          vel += node.vel();
        }
      }
    }
  }

  if(cnt > 0) vel /= cnt;
  else { 
    /* 
     * if the particle doesn't cover any nodes, we will set the velocity 
     * to the nodes the most cloest to particle. 
     */
    std::vector<double> 
      distance_t(_used_thread_num, std::numeric_limits<double>::max());
    std::vector<int> node_id_t(_used_thread_num, 0);

    int used_thread_num = _used_thread_num;
    if(_setting->sim_gridflow_enabled) {
      if(nodes_sum <= 400) used_thread_num = 1;
    }

    if(_setting->sim_gridflow_enabled) {
      #pragma omp parallel for num_threads(used_thread_num)
      for(int i=0; i<nodes_sum; i++) {
        
        int tid      = omp_get_thread_num();
        int node_id  = block_group.get(FlowBlockNodeType::flow_node, i);

        const auto & node = _chip->flow().nodes()[node_id];

        double dist = boost::geometry::distance(
          _particle->coord(), node.coord()
        );

        if(distance_t[tid] > dist) {
          distance_t[tid] = dist;
          node_id_t[tid]  = node_id;
        }
      }
    }
    else {
      #pragma omp parallel for
      for(size_t i = 0; i<_chip->flow().nodes().size(); i++) {
        int tid = omp_get_thread_num();

        const auto & node = _chip->flow().nodes()[i];
        double dist = boost::geometry::distance(
            _particle->coord(), node.coord());
        if(distance_t[tid] > dist) {
          distance_t[tid] = dist;
          node_id_t[tid]  = i;
        }
      }
    }

    // Merge results
    int node_id = node_id_t[0];
    double distance = distance_t[0];
    for(int i=0; i<used_thread_num; i++) {
      if(distance_t[i] < distance) {
        distance = distance_t[i];
        node_id  = node_id_t[i];
      }
    }

    if(covered_nodes_id) covered_nodes_id->push_back(node_id);
    vel = _chip->flow().node(node_id).vel();
  }

  return vel;
}

point 
Simulate::_apply_velocity(
  Particle * particle,
  const velocity & vel
)
{
  point coord = particle->coord();
  point new_coord(0, 0);

  new_coord = coord + vel * _setting->time_resolution;

  particle->update_coord(new_coord);
  return new_coord;
}

point 
Simulate::_apply_wall_effect(Particle *particle)
{
  point new_coord = particle->coord();
  bool collision = false;

  if(_setting->disable_wall_effect) {
    return new_coord;
  }

  if(_setting->adv_wall_effect_enabled) {
    collision = _apply_wall_effect_low_adv(particle, new_coord);
  }
  else {
    std::pair<point, point> two_points;
    double portion;
    collision = _apply_wall_effect_low(particle, two_points, portion);
    if(collision) {
      new_coord = _cal_wall_effect_position_based_on_A_B_porition(
        particle,
        two_points.first, two_points.second,
        portion
      );
    }
  }

  if(collision) {
    particle->update_coord(new_coord);

    Logger::add_record("Collision", to_string(particle->diameter()), 1);
  }

  return particle->coord();
}

point 
Simulate::_cal_wall_effect_position_based_on_A_B_porition(
  Particle *particle,
  point A,
  point B,
  double portion
)
{
  point middle       = (A + B) * 0.5;
  point direction    = particle->coord() - middle;
  point dist         = direction * (portion * _particle->radius());
  point new_position = particle->coord() + dist;

  return new_position;
}

bool
Simulate::_apply_wall_effect_low(
  /* input */
  Particle *particle,
  /* output */
  std::pair<point, point> &two_points,
  double &portion
)
{
  std::set<int> obstacles_cand;
  bool          collision = false;

  if(_setting->sim_gridflow_enabled) {
    flow_block_group  group_blocks;
    _flowgrids->get_adjacent_blocks(_particle, group_blocks);

    for(int i=0; i<group_blocks.size(FlowBlockNodeType::obstacle); i++) {
      obstacles_cand.insert(group_blocks.get(FlowBlockNodeType::obstacle, i));
    }
  }
  else {
    for(size_t i=0; i<_chip->design().obstacles().size(); i++) {
      obstacles_cand.insert(i);
    }
  }

  for(const auto & id : obstacles_cand) {

    const auto & obst = _chip->design().obstacles()[id];
    std::vector<point> res;
    if(particle->overlap_obstacle(obst, res, _debug_wall_effect)) {

      auto A = res.front();
      auto B = res.back();
      portion = (double) res.size() / 201.0;
      two_points = {A, B};
      collision = true;

      if(_debug_wall_effect) {
        point new_coord = _cal_wall_effect_position_based_on_A_B_porition(
          particle,
          A, B,
          portion
        );

        point middle       = (A + B) * 0.5;
        point direction    = particle->coord() - middle;

        std::cout << "!!!!" << portion << std::endl;
        std::cout << "Cov front = " << A << 
          " back = " << B << std::endl;
        std::cout << "Middle = " << middle << std::endl;
        std::cout << "Direction = " << direction << std::endl;

        std::cout << particle->coord() << " ===> " 
          << new_coord << std::endl;
        int a;
        std::cin >> a;
      }
    }
  }

  return collision;
}

bool 
Simulate::_apply_wall_effect_low_general_we(
  Particle *particle,
  point & new_coord
)
{
  point old_pos = particle->coord();
  bool debug = _is_debug_step();
  new_coord = old_pos;

  flow_block_group  group_blocks;

  GeneralWallEffect
    wall_effect(
      *particle, 
      *_flowgrids, 
      _chip->flow(), 
      group_blocks,
      _used_thread_num,
      particle->diameter() * _setting->collision_distance_ratio_threshold
    );

  bool ok = wall_effect.run();
  if(!ok) {
    printf("\n");
    LOG(WARNING) << "Timestep: " << _current_sim_step << " "
      << wall_effect.get_status_string() << '\n';
  }

  if(debug) {

    const auto & group = wall_effect.group();
    int nodes_sum = group.size(FlowBlockNodeType::obstacle_node);

    std::cout << "Calculate wall effect" << std::endl;
    std::cout << "\tobstacles nodes = " << nodes_sum << std::endl;
    std::cout << "\t# of on surface = " << wall_effect.num_on_surface_nodes() << std::endl;
    std::cout << "\t# of in particle = " << wall_effect.num_in_particle_nodes() << std::endl;
    std::cout << "\tdirection = " << wall_effect.adjsut_direction() << std::endl;
    std::cout << "\tdist_to_adjust = " << wall_effect.distance_adjustment() << std::endl;
    std::cout << "\tCollision = " << (wall_effect.has_collision()? "TRUE":"FALSE") 
      << std::endl;
    std::cout << "\told position = " << particle->coord() << std::endl;
    std::cout << "\tnew position = " << wall_effect.new_coord() << std::endl;

    std::filesystem::path path = _get_debug_file_path(
      ("debug_" + std::to_string(_current_sim_step) + "_"
       + to_string(_particle->diameter()) + "_obstacle.log")
    );
    _dump_debug_file(particle, false, true, path);

    std::filesystem::path path2 = _get_debug_file_path(
      ("debug_" + std::to_string(_current_sim_step) + "_" 
      + to_string(_particle->diameter()) + "_obstacle_adjustment.log")
    );

    std::ofstream out(path2);

    point direction = wall_effect.adjsut_direction();
    point adjustment = direction * wall_effect.distance_adjustment();

    out << adjustment.x() << ' '
        << adjustment.y() << '\n';

    out << wall_effect.new_coord().x() << ' '
        << wall_effect.new_coord().y() << '\n';

    out.close();
  }

  if(_is_debug_segment_step()) {
    const auto & group = wall_effect.group();
    int nodes_sum = group.size(FlowBlockNodeType::obstacle_node);
    for(int i=0; i<nodes_sum; i++) {
      int node_id = group.get(FlowBlockNodeType::obstacle_node, i);
      const auto & node = _chip->flow().nodes()[node_id];
      _debug_dump_segment->add_obstacle_nodes(node.coord());
    }
  }

  if(wall_effect.has_collision()) {
    new_coord = wall_effect.new_coord();

    if(_setting->dump_debug_collision_timestep) {
      std::filesystem::path output_path = _get_debug_file_path(
        _setting->chip_name + "_collision_timestep_" + 
        to_string(particle->diameter())
      );
      FILE *fp = fopen(output_path.c_str(), "a");
      fprintf(fp, "%d %lf\n", 
        _current_sim_step,
        wall_effect.distance_adjustment()
      );
      fclose(fp);
    }
  }

  return wall_effect.has_collision();
}

bool
Simulate::_apply_wall_effect_low_adv(
  Particle *particle,
  point & adv_new_position
)
{
  ASSERT(_setting->sim_gridflow_enabled == true);

  return _apply_wall_effect_low_general_we(particle, adv_new_position);

  flow_block_group group_blocks;
  int num_obstacle_nodes = 0;
  std::vector<point> points_on_circumference;
  std::vector<point> points_in_particle;
  std::vector<point> candidate_points;
  point adv_A(0, 0), adv_B(0, 0);
  bool adv_collision = false;
  bool adv_collision_with_candidate = false;;
  bool adv_collsion_by_candidate = false;

  double collision_threshold 
    = _setting->adv_wall_effect_collision_threshold;
  double candidate_points_threshold 
    = _setting->adv_wall_effect_candidate_threshold;
  double extra_shift 
    = _setting->adv_wall_effect_extra_shift;

  _flowgrids->get_adjacent_blocks(_particle, group_blocks);
  num_obstacle_nodes = group_blocks.size(FlowBlockNodeType::obstacle_node);

  for(int i=0; i<num_obstacle_nodes; i++) {
    int node_id = group_blocks.get(FlowBlockNodeType::obstacle_node, i);
    const auto & node = _chip->flow().node(node_id);

    double dist = boost::geometry::distance(_particle->coord(), node.coord());
    dist = (dist - particle->radius());

    if(fabs(dist) <= collision_threshold) {
      points_on_circumference.push_back(node.coord());
    }

    else if(dist > 0 and dist < candidate_points_threshold) {
      candidate_points.push_back(node.coord());
    }

    else if(dist < 0) {
      points_in_particle.push_back(node.coord());
    }
  }

  // Choose A B point
  if(points_on_circumference.size() == 1) {
    adv_B = adv_A = points_on_circumference.front();
    adv_collision = true;
  }
  else if(points_on_circumference.size() == 2) {
    adv_A = points_on_circumference.front();
    adv_B = points_on_circumference.back();
    adv_collision = true;
  }
  else if(points_on_circumference.size() > 3) {
    LOG(WARNING) << "Exceed two points on circumference.\n";

    adv_collision = true;
    std::sort(points_on_circumference.begin(), points_on_circumference.end(),
      [&](const point & a, const point &b)
      {
        point vec_a = a - _particle->coord();
        point vec_b = b - _particle->coord();

        double degree_a = std::atan2(vec_a.y(), vec_a.x());
        double degree_b = std::atan2(vec_b.y(), vec_b.x());

        return degree_a < degree_b;
      }
    );
    adv_A = points_on_circumference.front();
    adv_B = points_on_circumference.back();
  }
  else {
    if(points_in_particle.size()) {
      LOG(WARNING) << "Time resolution is too large."
        " Cannot detect wall during simulation. "
        "current step = \n" << _current_sim_step;
    }

    // Choose A B point from candidate_points
    if(candidate_points.size() > 1) {
      std::sort(candidate_points.begin(),candidate_points.end(),
        [&](const point & a, const point &b)
        {
          double dista = boost::geometry::distance(_particle->coord(), a);
          double distb = boost::geometry::distance(_particle->coord(), b);
          return dista < distb;
        }
      );
      point middle = 
        (candidate_points.front() + candidate_points.back())*0.5;

      double dist = boost::geometry::distance(_particle->coord(), middle);
      dist = (dist - particle->radius());

      if(fabs(dist) <= collision_threshold) {
        adv_collision = true;
        adv_collision_with_candidate = true;
        adv_A = adv_B = middle;
      }

      adv_collsion_by_candidate = true;
    }
  }

// Calculate new position
  if(adv_collision) {
    point middle = (adv_A + adv_B) * 0.5;
    point dir = _particle->coord() - middle;

    double len = boost::geometry::distance(dir, point(0,0));

    double dist = (_particle->radius() - len) + extra_shift;

    if(dist < 0) adv_collision = false;
    else adv_new_position = _particle->coord() + dir * (dist/len);
  }

  if(!adv_collision) adv_new_position = _particle->coord();

  bool consistency_check = false;
  consistency_check = 
    _setting->debug_adv_wall_effect_consistency | _debug_wall_effect;

// Consistency check
  if(consistency_check) {
    point error;
    bool adv_large_shift = false;
    Particle tmp_particle = *particle;
    Particle adv_particle = *particle;
    Particle golden_particle = *particle;

    adv_particle.update_coord(adv_new_position);

  // Calculate golden position
    std::pair<point, point> golden_two_points;
    point golden_A, golden_B;
    point golden_new_position;
    bool golden_collision = false;
    double golden_portion = 0;

    golden_collision = _apply_wall_effect_low(particle, 
      golden_two_points, golden_portion
    );

    golden_A = golden_two_points.first;
    golden_B = golden_two_points.second;

    golden_new_position = golden_collision? 
      _cal_wall_effect_position_based_on_A_B_porition(
        particle,
        golden_A, golden_B,
        golden_portion
      )
      :
      particle->coord()
    ;
    golden_particle.update_coord(golden_new_position);

    error = adv_new_position - golden_new_position;

    if(adv_collision)
      adv_large_shift = _is_large_shift(adv_new_position, particle->coord());

  // Check new position is collided.
    bool adv_new_pos_collision = false;
    if(adv_collision) {
      tmp_particle.update_coord(adv_new_position);
      double tmp_porition = 0;
      std::pair<point, point> tmp_two_points;
      adv_new_pos_collision = _apply_wall_effect_low(&tmp_particle, 
          tmp_two_points, tmp_porition);
    }

  // should we dump debug file
    bool debug = 
  (fabs(error.x()) >= 0.2 or fabs(error.y()) >= 0.2) 
        or adv_new_pos_collision
  or _debug_wall_effect
  or adv_large_shift
        or adv_collision_with_candidate 
  //or adv_A != adv_B
       //or golden_collision ^ adv_collision
    ;

    if(debug) 
    {
      std::cout << "\n";
      std::cout << "Particle coord = " << _particle->coord() << std::endl;
      std::cout << "Adv large shift = " << adv_large_shift << std::endl;
      std::cout << "Adv New position collision = " << adv_new_pos_collision 
  << std::endl;
      std::cout << "Adv collsion by candidate = " << 
        adv_collsion_by_candidate << std::endl;
      std::cout << "Adv points_on_circumference.size() = " <<
        points_on_circumference.size() << std::endl;
      std::cout << "Adv candidate_points.size() = " 
        << candidate_points.size() << std::endl;
      std::cout << "Adv collision with candidate points: " << 
        adv_collision_with_candidate << std::endl;

      std::cout << "Golden:\n";
      std::cout << "\tCollision : " << golden_collision << std::endl;
      std::cout << "\tA         : " << golden_A << std::endl;
      std::cout << "\tB         : " << golden_B << std::endl;
      std::cout << "\tNew coord : " << golden_particle.coord() << std::endl;
      std::cout << "Advance collision:\n";
      std::cout << "\tCollision : " << adv_collision << std::endl;
      std::cout << "\tA         : " << adv_A << std::endl;
      std::cout << "\tB         : " << adv_B << std::endl;
      std::cout << "\tNew coord : " << adv_new_position << std::endl;

      std::cout << "===================================\n";
      std::cout << "Error: " <<  error << std::endl << std::endl;

      point middle = (adv_A + adv_B) * 0.5;
      std::cout << "dist(Particle,middle) = " 
        << boost::geometry::distance(middle, _particle->coord()) << std::endl;

      for(auto point : candidate_points) {
        double dista = boost::geometry::distance(_particle->coord(), point);
        std::cout << dista << std::endl;
      }

      if(adv_collision) 
      {
        if(adv_collsion_by_candidate == false) {
          point middle = (adv_A + adv_B) * 0.5;
          double dist = boost::geometry::distance(middle, _particle->coord());
          std::cout << "dist(Particle, middle) = " << dist << std::endl;
          std::cout << "Move distacnce = " 
            << _particle->radius() - dist << std::endl;
        }
        else {
          point cand_middle = 
            (candidate_points.front() + candidate_points.back()) *0.5;
          std::cout << "dist(Particle,cand_middle) = " 
            << boost::geometry::distance(cand_middle, _particle->coord()) << std::endl;
        }
      }

      // dump debug information
      std::ofstream out("debug_advance_collision.txt");
      particle->dump(out);
      adv_particle.dump(out);
      golden_particle.dump(out);
      group_blocks.dump_nodes(out, FlowBlockNodeType::obstacle_node, 
        [this](std::ostream& os, int id) {
          os << _chip->flow().node(id).x() << ' '
             << _chip->flow().node(id).y();
        }
      );
      out << golden_collision << std::endl;
      out << golden_A.x() << " " << golden_A.y() << std::endl;
      out << golden_B.x() << " " << golden_B.y() << std::endl;
      out << adv_collision << std::endl;
      out << adv_A.x() << " " << adv_A.y() << std::endl;
      out << adv_B.x() << " " << adv_B.y() << std::endl;
      out.close();

      //for(auto [dist, id] : nodes_distance) {
  //std::cout << dist << ", " << id << "\n";
      //}
      int a;
      std::cin >> a;
    }
  }

  return adv_collision;
}

//double 
//Simulate::_cal_adv_wall_effect_portion(
  //Particle *particle,
  //point A,
  //point B
//)
//{
  //if(A == B) {
    //double distance = boost::geometry::distance(particle->coord(), A);
    //double diff = distance - particle->radius();
    //if(diff < 0) {
      //return std::fabs(diff) / (particle->radius() * std::fabs(distance));
    //}
    //return 1.0/201;
  //}
  //else {
    ////TODO:
    //Logger::add_record("Adv collision", "Unable to calculate portion", 1);
    //return 1.0/201;
  //}
//}

bool 
Simulate::_is_large_shift(
  point new_coord, 
  point old_coord
)
{
  point diff = new_coord - old_coord;
  return (std::fabs(diff.x()) >= 0.3) || (std::fabs(diff.y()) >= 0.3);
}

bool 
Simulate::_is_debug_step()
{
  if(_debug_step.size() == 0) return false;
  else return _debug_step.find(_current_sim_step) != _debug_step.end();
}

void 
Simulate::_dump_trajectory(
  const std::vector<point> &trajectory,
  const std::filesystem::path &output_path
)
{
  LOG(INFO) << "Dump result to " << output_path << "\n";
  FILE *fp = fopen(output_path.c_str(), "w");
  fprintf(fp, "%ld %f\n", trajectory.size(), _particle->diameter());
  for(const auto & point : trajectory) {
    fprintf(fp, "%lf %lf\n", point.x(), point.y());
  }
  fclose(fp);
}

std::filesystem::path 
Simulate::_get_debug_file_path(const std::string &filename)
{
  std::filesystem::path debug_folder = _setting->debug_output_folder;
  return debug_folder / filename;
}

bool 
Simulate::_is_debug_segment_step() const
{
  return _debug_dump_segment != nullptr 
  && (_setting->debug_range_start <= _current_sim_step)
  && (_current_sim_step <= _setting->debug_range_stop);
}

void
Simulate::_dump_debug_file(
  Particle *particle,
  bool include_vel, 
  bool include_obstacle,
  const std::filesystem::path &output_path
) const
{
  std::ofstream out(output_path);

  out << std::fixed << std::setprecision(6);

  out << particle->coord().x() << ' ' 
      << particle->coord().y() << '\n';
  out << particle->diameter() << std::endl;

  flow_block_group group;
  _flowgrids->get_adjacent_blocks(_particle, group);

  if(include_vel) {
    group.dump_nodes(out, FlowBlockNodeType::flow_node, 
      [this](std::ostream &os, int id) {
        const auto & node = _chip->flow().node(id);
        os << node.x() << ' '
           << node.y() << ' '
           << _particle->cover(node.coord());
      }
    );
  }
  else out << "0\n";

  if(include_obstacle) {
    group.dump_nodes(out, FlowBlockNodeType::obstacle_node, 
      [this, particle](std::ostream &os, int id) {
        const auto & node = _chip->flow().node(id);
        bool on_surface, in_particle;

        double distance = boost::geometry::distance(
          particle->coord(), node.coord()
        );
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

        os << node.x() << ' '
           << node.y() << ' '
           << on_surface << ' '
           << in_particle;
      }
    );
  }
  else out << "0\n";

  out.close();
}

}
