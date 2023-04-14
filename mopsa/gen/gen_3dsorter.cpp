#include <mopsa/logger/logger.hpp>
#include <mopsa/gen/gen_3dsorter.hpp>
#include <mopsa/sim/sim.hpp>

#include <random>

using namespace mopsa::Gen3DSorter;
using namespace mopsa;

//  Class : GenSetting
//-----------------------------------------------------------------------------
GenSetting::GenSetting()
{
  _register_variable(
    "output_folder",
    SettingVarType::String,
    gen_path
  );

  _register_variable(
    "target_particle_dims",
    SettingVarType::FloatVec, 
    target_particle_dims
  );

  _register_variable(
    "ga_num_elite",
    SettingVarType::Int, 
    ga_num_elite 
  );

  _register_variable(
    "ga_initial_population_size",
    SettingVarType::Int, 
    ga_init_population_size
  );

  _register_variable(
    "ga_population_size",
    SettingVarType::Int, 
    ga_population_size
  );

  _register_variable(
    "ga_num_iterations",
    SettingVarType::Int, 
    ga_num_iterations
  );
  
  _register_variable(
    "ga_keep_elite_to_next_generation",
    SettingVarType::Boolean,
    ga_keep_elite_to_next_generation
  );

  _register_variable(
    "recipe_radius_range_min_in_um",
    SettingVarType::Float,
    recipe_radius_range_min_in_um
  );

  _register_variable(
    "recipe_radius_range_max_in_um",
    SettingVarType::Float,
    recipe_radius_range_max_in_um
  );

  _register_variable(
    "recipe_theta_range_min_in_degree",
    SettingVarType::Float,
    recipe_theta_range_min_in_degree
  );

  _register_variable(
    "recipe_theta_range_max_in_degree",
    SettingVarType::Float,
    recipe_theta_range_max_in_degree
  );

  _register_variable(
    "recipe_phi_range_min_in_degree",
    SettingVarType::Float,
    recipe_phi_range_min_in_degree
  );

  _register_variable(
    "recipe_phi_range_max_in_degree",
    SettingVarType::Float,
    recipe_phi_range_max_in_degree
  );

  _register_variable(
    "sim_samples_count",
    SettingVarType::Int,
    monte_carlo_samples_num
  );

  _register_variable(
    "sim_3sigma_radius",
    SettingVarType::Float,
    monte_carlo_3sigma_radius
  );
}

void
GenSetting::dump(std::ostream& os) 
{
  os << "*********************************\n";
  os << "   Design Automation Setting     \n";
  os << "*********************************\n";
  os << "output_folder: " << gen_path << '\n';
  os << "target_particle_dims: "; 
    _dump_vec(os, target_particle_dims); os << '\n';

  os << "ga_num_elite: " << ga_num_elite << '\n';
  os << "ga_init_population_size: " << ga_init_population_size << '\n';
  os << "ga_population_size: " << ga_population_size << '\n';
  os << "ga_num_iterations: " << ga_num_iterations << '\n';
  os << "recipe_radius_range: [" << recipe_radius_range_min_in_um << 
    ", " << recipe_radius_range_max_in_um << "]\n";
  os << "recipe_theta_range: [" << recipe_theta_range_min_in_degree << 
    ", " << recipe_theta_range_max_in_degree << "]\n";
  os << "recipe_phi_range: [" << recipe_phi_range_min_in_degree << 
    ", " << recipe_phi_range_max_in_degree << "]\n";
  os << "sim_samples_count: " << monte_carlo_samples_num << '\n';
  os << "sim_3sigma_radius: " << monte_carlo_3sigma_radius << '\n';
}

//  Class : GenManager
//-----------------------------------------------------------------------------
Recipe
GenManager::random_recipe(GenSetting *setting)
{
  Recipe recipe;

  recipe.radius = rand_range(
    setting->recipe_radius_range_min_in_um, 
    setting->recipe_radius_range_max_in_um
  );

  recipe.y_init = rand_range(
    0,
    recipe.radius
  );

  recipe.y_init_shift = rand_range(
    0,
    recipe.radius
  );

  recipe.y_gap = rand_range(
    setting->recipe_radius_range_min_in_um, 
    setting->recipe_radius_range_max_in_um
  );

  recipe.x_init = rand_range(
    0,
    recipe.radius
  );

  recipe.x_gap = rand_range(
    setting->recipe_radius_range_min_in_um, 
    setting->recipe_radius_range_max_in_um
  );

  recipe.theta = rand_range(
    setting->recipe_theta_range_min_in_degree, 
    setting->recipe_theta_range_max_in_degree
  );

  recipe.phi = rand_range(
    setting->recipe_phi_range_min_in_degree, 
    setting->recipe_phi_range_max_in_degree
  );

  recipe.fixed_side = setting->recipe_fixed_side;

  return recipe;
}

bool 
GenManager::launch_matlab(
  GenSetting* setting, 
  const GenDesign &design
)
{
  const auto & matlab_setting_path = 
    std::filesystem::absolute(design.get_design_folder()/ "matlab_setting.m");

  FILE *fp = fopen(matlab_setting_path.c_str(), "w");
  fprintf(fp, "s_recipe_path = \'%s\';\n", 
    std::filesystem::absolute(design.get_recipe_path()).c_str()
  );
  fprintf(fp, "s_model_name = \'%s\';\n", design.get_design_name().c_str());
  fprintf(fp, "s_model_path = \'%s\';\n", 
    std::filesystem::absolute(design.get_design_folder()).c_str()
  );
  fclose(fp);

  const auto& cur_path = std::filesystem::current_path();

  std::string launch_matlab = 
    setting->run_matlab_script_path
    + " " 
    + std::string(matlab_setting_path)
    + " > "
    + std::string(
        std::filesystem::absolute(design.get_design_folder() / "comsol.log") 
      )
  ;

  std::string total_command = 
    launch_matlab + ";"
    + "cd " + std::string(cur_path) + ";"
  ;

  LOG(INFO) << "Build comsol " << design.get_design_name() 
    << ": "<< total_command << '\n';
  std::cout << total_command << std::endl;

  return std::system(total_command.c_str()) == 0;
}

bool 
GenManager::launch_simulation(
  GenSetting* setting, 
  const GenDesign &design
)
{
  // sim_setting
  const auto & sim_setting = 
    GenManager::get_simulation_setting(setting, design);

  const auto &sim_setting_path = 
    GenManager::get_sim_setting_path(design.get_design_folder());

  GenManager::dump_sim_setting(sim_setting, sim_setting_path);

  const auto& cur_path = std::filesystem::current_path();

  std::string cd_cmd = "cd " + std::string(design.get_design_folder()) + "; ";
  std::string rm_cmd = "rm -rf " + std::string(
      std::filesystem::absolute(design.get_design_folder()/"sim.log")) + "; ";
  std::string launch_cmd = 
    setting->run_evaluate_script_path
    + " " 
    + std::string(sim_setting_path)
    + " > "
    + std::filesystem::absolute(design.get_design_folder() / "sim.log").string()
    + "; "
  ;

  std::string total_command = 
      cd_cmd 
    + rm_cmd 
    + launch_cmd
    + " cd " + std::string(cur_path) + ";"
  ;

  LOG(INFO) << "Evaluate design " << design.get_design_name() 
    << ": "<< total_command << '\n';
  return std::system(total_command.c_str()) == 0;
}

float
GenManager::rand_range(float begin, float end) 
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dist(begin, end);

  return dist(gen);
};

std::filesystem::path 
GenManager::get_gen_folder(GenSetting *setting, int gid)
{
  return std::filesystem::path(setting->gen_path) / ("G" + std::to_string(gid));
}

std::string
GenManager::get_design_name(GenSetting* setting, int gen_id, int design_id)
{
  return "G" + std::to_string(gen_id) + "_" + std::to_string(design_id);
}

std::filesystem::path
GenManager::get_simulation_points_path(GenSetting *setting)
{
  return std::filesystem::absolute(
    std::filesystem::path(setting->gen_path) / "simulation_startpoints.txt"
  );
}

std::filesystem::path
GenManager::get_recipe_path(const std::filesystem::path & design_folder)
{
  return std::filesystem::absolute(
    design_folder / "feature_recipe.txt"
  );
}

std::filesystem::path
GenManager::get_velcoity_path(const std::filesystem::path & design_folder)
{
  return std::filesystem::absolute(
    design_folder / "vel_size_3.csv"
  );
}

std::filesystem::path
GenManager::get_score_path(const std::filesystem::path & design_folder)
{
  return std::filesystem::absolute(
    design_folder / "sim_output/score.txt"
  );
}

std::filesystem::path
GenManager::get_sim_setting_path(const std::filesystem::path & design_folder)
{
  return std::filesystem::absolute(
    design_folder / "sim_setting.m"
  );
}

float
GenManager::get_simulation_startpoint_z(GenSetting *setting)
{
  return setting->model_container_height/2 + setting->model_input_region_height;
}

float
GenManager::get_simulation_boundary_z(GenSetting *setting)
{
  return -setting->model_container_height/2 
    - setting->model_output_region_height + 5;
}

SimSetting3d
GenManager::get_simulation_setting(
  GenSetting *setting, 
  const GenDesign &design
)
{
  SimSetting3d sim_setting;

  sim_setting.chip_name = design.get_design_name();
  sim_setting.mesh_nodes_path = design.get_velcoity_nodes_path();
  sim_setting.output_folder = std::filesystem::absolute(
    design.get_design_folder() / "sim_output");

  sim_setting.dPs = setting->target_particle_dims;

  //sim_setting.boundary_max_timestep = 1000000;
  sim_setting.boundary_max_timestep = 500000;
  sim_setting.boundary_less_z = GenManager::get_simulation_boundary_z(setting);
  sim_setting.time_resolution = 1;

  sim_setting.init_position_path = 
    GenManager::get_simulation_points_path(setting);

  sim_setting.num_threads = 16;

  sim_setting.dump_detail_trajectory = false;
  return sim_setting;
}

void
GenManager::dump_sim_setting(
  const SimSetting3d & sim_setting, 
  const std::filesystem::path &path
)
{
  FILE *fp = fopen(path.c_str(), "w");
  fprintf(fp, "chip_name = \"%s\";\n", 
    sim_setting.chip_name.c_str()
  );
  fprintf(fp, "mesh_nodes_path = \"%s\";\n", 
    sim_setting.mesh_nodes_path.c_str()
  );
  fprintf(fp, "output_folder = \"%s\";\n", 
    sim_setting.output_folder.c_str()
  );

  fprintf(fp, "dPs = [%lf", sim_setting.dPs.front());
  for(size_t i=1; i<sim_setting.dPs.size(); i++) {
    fprintf(fp, ", %lf", sim_setting.dPs[i]);
  }
  fprintf(fp, "];\n");

  fprintf(fp, "time_resolution = %lf;\n",
    sim_setting.time_resolution
  );
  fprintf(fp, "boundary_max_timestep = %d;\n",
    sim_setting.boundary_max_timestep
  );
  fprintf(fp, "boundary_less_z = %lf;\n", 
    sim_setting.boundary_less_z
  );
  fprintf(fp, "init_position_path = \"%s\";\n", 
    sim_setting.init_position_path.c_str()
  );
  fprintf(fp, "num_threads = %d;\n",
    sim_setting.num_threads
  );
  fprintf(fp, "dump_detail_trajectory = %d;\n",
    sim_setting.dump_detail_trajectory
  );

  fclose(fp);
}

//  Class : GeneticOpt
//-----------------------------------------------------------------------------

GeneticOpt::GeneticOpt(GenSetting *setting)
  :
  _setting(setting)
{
  _current_generation_id = -1;
}

bool
GeneticOpt::run()
{
  if(!std::filesystem::exists(_setting->gen_path)) {
    std::filesystem::create_directory(_setting->gen_path);
  }

  if(_current_generation_id == -1) {
    _gen_monte_carlo_samples();

    _gen_first_generation();

    _current_generation_id = 0;
  }
  else {
    // make sure the current generation is evaluated.
    update_generation(_current_generation_id);
  }

  for(int i=_current_generation_id; i<_setting->ga_num_iterations; i++) {

    LOG(INFO) << "Generate generation: " << i + 1 << '\n';
    _gen_next_generation(i);

    _current_generation_id += 1;
  }

  return true;
}

bool
GeneticOpt::restore_session(
  const std::filesystem::path &path
)
{
  std::filesystem::directory_iterator dirs{path};

  for(const auto & entry: std::filesystem::directory_iterator{path}) {

    if(int gid = -1; 
      entry.is_directory() && _is_generation_folder(entry.path(), gid)) 
    {
      _restore_generation_folder(entry.path(), gid);
    }
  }

  const auto & sim_startpoints_path =
    GenManager::get_simulation_points_path(_setting);
  if(std::filesystem::exists(sim_startpoints_path)) 
  {
    _restore_monte_carlo_samples(sim_startpoints_path);
  }
  else {
    _gen_monte_carlo_samples();
  }

  LOG(INFO) << "Loaded " << _current_generation_id + 1 << " generations.\n";

  return true;
}

void 
GeneticOpt::_gen_monte_carlo_samples()
{
  float fixed_z = GenManager::get_simulation_startpoint_z(_setting);

  std::default_random_engine generator;

  std::normal_distribution<double> distribution_x(
    _setting->monte_carlo_center.x(), _setting->monte_carlo_3sigma_radius/3.0
  );

  std::normal_distribution<double> distribution_y(
    _setting->monte_carlo_center.y(), _setting->monte_carlo_3sigma_radius/3.0
  );

  for(int i=0; i<_setting->monte_carlo_samples_num; i++) {
    float x = distribution_x(generator);
    float y = distribution_y(generator);
    _simulating_points.emplace_back(x, y, fixed_z);
  }

  const auto &sim_startpoints_path = 
    GenManager::get_simulation_points_path(_setting);

  FILE *fp = fopen(sim_startpoints_path.c_str(), "w");
  fprintf(fp, "%ld\n", _simulating_points.size());
  for(const auto &point : _simulating_points) {
    fprintf(fp, "%f %f %f\n", point.x(), point.y(), point.z());
  }
  fclose(fp);
}

void 
GeneticOpt::_gen_first_generation()
{
  constexpr int gen_id = 0;

  const auto & cur_generation_folder = 
    GenManager::get_gen_folder(_setting, gen_id);

  std::filesystem::create_directory(cur_generation_folder);

  for(int i=0; i<_setting->ga_init_population_size; i++) {

    const auto &design_name = GenManager::get_design_name(_setting, gen_id, i);
    const auto &design_folder = cur_generation_folder / design_name;
    const auto &recipe = GenManager::random_recipe(_setting);
    const auto &recipe_path = GenManager::get_recipe_path(design_folder);

    _get_design(gen_id, i).init(
      gen_id,
      i,
      design_name,
      design_folder,
      recipe_path,
      recipe
    );
  }

  _build_comsol_model(gen_id);

  _evaluate_design(gen_id);
}

void
GeneticOpt::_gen_next_generation(int gid)
{
  int next_gid = gid + 1;
  const auto & next_generation_folder = 
    GenManager::get_gen_folder(_setting, next_gid);
  std::filesystem::create_directory(next_generation_folder);

  const auto &cur_designs = get_sorted_design(gid);
  int next_design_id = 0;
  std::vector<GenDesign*> parents;
  int num_elite = std::min(_setting->ga_num_elite, (int)cur_designs.size());
  for(int i = 0; i < num_elite; i++) {
    parents.push_back(cur_designs[i]);

    if(!_setting->ga_keep_elite_to_next_generation) continue;

    const auto &design_name = 
      GenManager::get_design_name(_setting, next_gid, next_design_id);
    const auto &design_folder = next_generation_folder / design_name;
    const auto &recipe = cur_designs[i]->get_recipe();
    const auto &recipe_path = GenManager::get_recipe_path(design_folder);

    auto &new_design = _get_design(next_gid, next_design_id);

    new_design.init(
      next_gid,
      next_design_id,
      design_name,
      design_folder,
      recipe_path,
      recipe
    );

    LOG(INFO) << "Copy design: " << cur_designs[i]->get_design_name() << " to "
      << new_design.get_design_name() << '\n';
    new_design.copy_result_from(*cur_designs[i]);

    const auto & sim_setting = 
      GenManager::get_simulation_setting(_setting, new_design);
    const auto &sim_setting_path = 
      GenManager::get_sim_setting_path(new_design.get_design_folder());
    GenManager::dump_sim_setting(sim_setting, sim_setting_path);

    next_design_id++;
  }

  std::vector<Recipe> from_crossover;
  std::vector<std::string> logs;
  while( 
    (next_design_id + (int)from_crossover.size())
    < _setting->ga_population_size
  ) 
  {
    auto [d1, d2] = _randomly_select_two_designs(parents);

    std::cout << d1->get_design_name() << "(" << d1->get_score() << ") ";
    std::cout << d2->get_design_name() << "(" << d2->get_score() << ")\n";

    const auto [r1, r2] = _random_merge_two_recipe_v2(
      d1->get_recipe(), 
      d2->get_recipe()
    );

    std::string log = "parent1: " + d1->get_design_name() 
      + ", " + std::to_string(d1->get_score()) + "\n"
      + "parent2: " + d2->get_design_name() + ", "
      + std::to_string(d2->get_score()) + "\n";

    from_crossover.push_back(r1);
    from_crossover.push_back(r2);

    logs.push_back(log);
    logs.push_back(log);
  }

  for(int i=0; 
    next_design_id < _setting->ga_population_size; 
    i++, next_design_id++
  ) 
  {
    Recipe recipe = from_crossover[i];
    std::string log;
    _recipe_mutation(recipe, log);

    std::string final_log = logs[i] + log;

    auto &new_design = _get_design(next_gid, next_design_id);
    const auto &design_name = 
      GenManager::get_design_name(_setting, next_gid, next_design_id);
    const auto &design_folder = next_generation_folder / design_name;
    const auto &recipe_path = GenManager::get_recipe_path(design_folder);

    new_design.init(
      next_gid,
      next_design_id,
      design_name,
      design_folder,
      recipe_path,
      recipe
    );

    const auto & evolution_log = design_folder / "evolution.log";
    FILE *fp = fopen(evolution_log.c_str(), "w");
    fprintf(fp, "%s\n", final_log.c_str());
    fclose(fp);
  }

  _build_comsol_model(next_gid);

  _evaluate_design(next_gid);
}

void
GeneticOpt::update_generation(int gid) 
{
  LOG(INFO) << "Update generation: " << gid << '\n';

  _build_comsol_model(gid);
  _evaluate_design(gid);
}

Recipe
GeneticOpt::_random_merge_two_recipe(const Recipe& r1, const Recipe &r2)
{
  Recipe recipe = r1;

  if(rand()%2) recipe.radius = r2.radius;
  if(rand()%2) recipe.theta = r2.theta;
  if(rand()%2) recipe.phi = r2.phi;

  if(rand()%2) recipe.y_init = r2.y_init;
  if(rand()%2) recipe.y_init_shift = r2.y_init_shift;
  if(rand()%2) recipe.y_gap = r2.y_gap;

  if(rand()%2) recipe.x_init = r2.x_init;

  if(rand()%2) recipe.fixed_side = r2.fixed_side;


  // _legalize_recipe(recipe);

  return recipe;
}

std::pair<Recipe, Recipe>
GeneticOpt::_random_merge_two_recipe_v2(const Recipe& r1, const Recipe &r2)
{
  Recipe recipe1;
  Recipe recipe2;

  if(rand()%2) {
    recipe1.radius = r1.radius;
    recipe2.radius = r2.radius;
  }
  else {
    recipe1.radius = r2.radius;
    recipe2.radius = r1.radius;
  }

  if(rand()%2) {
    recipe1.theta = r1.theta;
    recipe2.theta = r2.theta;
  }
  else {
    recipe1.theta = r2.theta;
    recipe2.theta = r1.theta;
  }

  if(rand()%2) {
    recipe1.phi = r1.phi;
    recipe2.phi = r2.phi;
  }
  else {
    recipe1.phi = r2.phi;
    recipe2.phi = r1.phi;
  }

  if(rand()%2) {
    recipe1.y_init = r1.y_init;
    recipe2.y_init = r2.y_init;
  }
  else {
    recipe1.y_init = r2.y_init;
    recipe2.y_init = r1.y_init;
  }

  if(rand()%2) {
    recipe1.y_init_shift = r1.y_init_shift;
    recipe2.y_init_shift = r2.y_init_shift;
  }
  else {
    recipe1.y_init_shift = r2.y_init_shift;
    recipe2.y_init_shift = r1.y_init_shift;
  }

  if(rand()%2) {
    recipe1.y_gap = r1.y_gap;
    recipe2.y_gap = r2.y_gap;
  }
  else {
    recipe1.y_gap = r2.y_gap;
    recipe2.y_gap = r1.y_gap;
  }

  if(rand()%2) {
    recipe1.x_init = r1.x_init;
    recipe2.x_init = r2.x_init;
  }
  else {
    recipe1.x_init = r2.x_init;
    recipe2.x_init = r1.x_init;
  }

  if(rand()%2) {
    recipe1.x_gap = r1.x_gap;
    recipe2.x_gap = r2.x_gap;
  }
  else {
    recipe1.x_gap = r2.x_gap;
    recipe2.x_gap = r1.x_gap;
  }

  if(rand()%2) {
    recipe1.fixed_side = r1.fixed_side;
    recipe2.fixed_side = r2.fixed_side;
  }
  else {
    recipe1.fixed_side = r2.fixed_side;
    recipe2.fixed_side = r1.fixed_side;
  }

  return {recipe1, recipe2};
}

void
GeneticOpt::_legalize_recipe(Recipe &recipe)
{
  if(recipe.y_init >= recipe.radius) {
    recipe.y_init = recipe.radius;
  }

  if(recipe.y_init >= recipe.radius) {
    recipe.y_init = recipe.radius;
  }

  if(recipe.x_init >= recipe.radius) {
    recipe.x_init = recipe.radius;
  }
}

void
GeneticOpt::_recipe_mutation(Recipe& r1, std::string &log)
{
  if(rand()%100 >= 5) {
    log = "No mutation\n";
    return;
  }

  Recipe r2 = GenManager::random_recipe(_setting);

  size_t num_mutation = rand() % 3 + 1;

  std::set<int> to_mutation;
  while(to_mutation.size() < num_mutation) {
    to_mutation.insert(rand()%8);
  }

  log = "Num_numtation: " + std::to_string(to_mutation.size()) + '\n';
  for(auto x : to_mutation) {
    log += std::to_string(x) + "; ";
  }
  log += '\n';

  if(to_mutation.find(0) != to_mutation.end()) {
    r1.radius = r2.radius;
  }

  if(to_mutation.find(1) != to_mutation.end()) {
    r1.theta = r2.theta;
  }

  if(to_mutation.find(2) != to_mutation.end()) {
    r1.phi = r2.phi;
  }

  if(to_mutation.find(3) != to_mutation.end()) {
    r1.y_init = r2.y_init;
  }

  if(to_mutation.find(4) != to_mutation.end()) {
    r1.y_init_shift = r2.y_init_shift;
  }

  if(to_mutation.find(5) != to_mutation.end()) {
    r1.y_gap = r2.y_gap;
  }

  if(to_mutation.find(6) != to_mutation.end()) {
    r1.x_init = r2.x_init;
  }

  if(to_mutation.find(7) != to_mutation.end()) {
    r1.x_gap = r2.x_gap;
  }

  // _legalize_recipe(recipe);
}

std::pair<GenDesign*, GenDesign*>
GeneticOpt::_randomly_select_two_designs(
  const std::vector<GenDesign*> &parents
)
{
  int id1 = (int)GenManager::rand_range(0, parents.size());
  int id2 = (int)GenManager::rand_range(0, parents.size());
  
  while(id2 == id1) {
    id2 = (int)GenManager::rand_range(0, parents.size());
  }

  return {
    parents[id1],
    parents[id2],
  };
}

void
GeneticOpt::_build_comsol_model(int gid) 
{
  LOG(INFO) << "Build comsol model for generation: " << gid << '\n';

  auto & cur_gen_designs = _designs[gid];

  for(size_t i=0; i<cur_gen_designs.size(); i++) {

    auto & design = _get_design(gid, i);
    
    if(design.get_status() != GenDesignStatus::BUILT_RECIPE) {
      LOG(WARNING) << "Skip " << design.get_design_name() << " in comsol. "
        << "Status = " << design.get_status() << "\n";
      continue;
    }

    // call matlab to build comosl model
    GenManager::launch_matlab(_setting, design);

    design.update_vel_path();

    LOG(INFO) << "Design: " << design.get_design_name() << ": " << 
      design.get_status() << '\n';
  }
}

void
GeneticOpt::_evaluate_design(int gid)
{
  LOG(INFO) << "Evaluation for generation: " << gid << '\n';

  auto & cur_gen_designs = _designs[gid];

  for(size_t i=0; i<cur_gen_designs.size(); i++) {

    auto & design = _get_design(gid, i);
    ASSERT(!design.is_uninit());

    if(design.get_status() != GenDesignStatus::BUILT_COMSOL) {
      LOG(WARNING) << "Skip " << design.get_design_name() << " in evaluation. "
        << "Status = " << design.get_status() << "\n";
      continue;
    }

    bool ok = GenManager::launch_simulation(_setting, design);

    if(ok) ok = design.update_score();

    if(ok) {
      ASSERT(design.get_status() == GenDesignStatus::SIMULATED);
      LOG(INFO) << "evalulated successfully, score: " 
        << design.get_score() << '\n';
    }
  }
}

bool 
GeneticOpt::_is_generation_folder(
  const std::filesystem::path &path, 
  int &gid
)
{
  const auto & name = std::string(path.stem());
  if(name[0] == 'G') {
    gid = 0;
    for(size_t i=1; i<name.size(); i++) {
      if(name[i] == '_') break;
      ASSERT(name[i] >= '0' and name[i] <= '9');
      gid = gid * 10 + name[i] - '0';
    }
    return true;
  }
  else {
    gid = -1;
    return false;
  }
}

bool 
GeneticOpt::_restore_generation_folder(
  const std::filesystem::path &path,
  int gen_id
)
{
  LOG(INFO) << "Restore generation_folder: " << path 
    << " gid = " << gen_id << '\n';

  const auto & gen_base_name = std::string(path.stem());

  for(const auto & entry: std::filesystem::directory_iterator{path}) {
    if(entry.is_directory()) {
      const auto &folder_name = std::string(entry.path().stem());
      if(folder_name.find(gen_base_name) == 0) {
        int design_id = 0;
        for(size_t i=gen_base_name.size()+1; i<folder_name.size(); i++) {
          ASSERT(folder_name[i] >= '0' and folder_name[i] <= '9');
          design_id = design_id * 10 + folder_name[i] - '0';
        }

        if(!_restore_design_folder(entry.path(), gen_id, design_id)) {
          return false;
        }
      }
    }
  }

  _current_generation_id = std::max(_current_generation_id, gen_id);

  return true;
}

bool 
GeneticOpt::_restore_design_folder(
  const std::filesystem::path &design_path, 
  int gen_id, 
  int design_id
)
{
  LOG(INFO) << "Restore design: " << design_path.stem()
    << " from " << design_path 
    << " generation_id " << gen_id 
    << " design_id " << design_id 
  << '\n';

  auto &design = _get_design(gen_id, design_id);

  bool ok = design.restore_session(
    gen_id,
    design_id,
    design_path.stem(),
    design_path
  );

  LOG(INFO) << "=> " << design.get_status() 
    << ' ' << design.get_score() << '\n';

  return ok;
}

bool 
GeneticOpt::_restore_monte_carlo_samples(const std::filesystem::path &path)
{
  LOG(INFO) << "Restore simulation startpoints: " << path << '\n';

  FILE *fp = fopen(path.c_str(), "r");
  int num = 0;
  ASSERT(fscanf(fp, "%d\n", &num) == 1);
  ASSERT(num == _setting->monte_carlo_samples_num);
  _simulating_points.resize(num);
  for(int i=0; i<num; i++) {
    float x, y, z;
    ASSERT( fscanf(fp, "%f %f %f\n", &x, &y, &z) == 3);
    _simulating_points[i] = {x, y, z};
  }
  fclose(fp);

  LOG(INFO) << num << " points are loaded\n";
  return true;
}

GenDesign&
GeneticOpt::_get_design(
  size_t gid, 
  size_t did
) 
{
  if(_designs.size() < gid + 1) {
    _designs.resize(gid + 1);
  }

  if(_designs[gid].size() < did + 1) {
    _designs[gid].resize(did + 1);
  }

  return _designs[gid][did];
}

std::vector<GenDesign*> 
GeneticOpt::get_sorted_design(int gid)
{
  auto & cur_gen_designs = _designs[gid];

  std::vector<GenDesign*> designs;

  for(size_t i=0; i<cur_gen_designs.size(); i++) {
    designs.push_back(&cur_gen_designs[i]);
  }

  std::sort(
    designs.begin(), designs.end(), 
    [](GenDesign *d1, GenDesign *d2) {
      return d1->get_score() > d2->get_score();
    }
  );

  return designs;
}
