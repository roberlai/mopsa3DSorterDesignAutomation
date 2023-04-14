#ifndef MOPSA_GEN_3DSORTER_H
#define MOPSA_GEN_3DSORTER_H

#include <mopsa/headerdef.hpp>
#include <mopsa/gen/gen_3design.hpp>
#include <mopsa/util/reader.hpp>
#include <mopsa/sim/sim.hpp>

#include <filesystem>

namespace mopsa {

namespace Gen3DSorter {

//  Class : GenSetting
//-----------------------------------------------------------------------------
class GenSetting : public SettingBase
{

public:

  GenSetting();

/*------ parameters ------*/
  const std::string run_matlab_script_path{
    "${MOPSA_ROOT}/scripts/runMopsaAuto3DSorter.sh"
  };
  const std::string run_evaluate_script_path{
    "${MOPSA_ROOT}/scripts/evaluate3D.sh"
  };
  std::string gen_path{"Gen3dSorter"};

  std::vector<double> target_particle_dims{5, 15};
  int ga_num_elite{100};
  int ga_init_population_size{500};
  int ga_population_size{1000};
  int ga_num_iterations{10};
  bool ga_keep_elite_to_next_generation{false};

  double recipe_radius_range_min_in_um{50};
  double recipe_radius_range_max_in_um{60};
  double recipe_theta_range_min_in_degree{30};
  double recipe_theta_range_max_in_degree{150};
  double recipe_phi_range_min_in_degree{0};
  double recipe_phi_range_max_in_degree{40};
  const int recipe_fixed_side{1}; // fixed at x by default. (change y and z)

  int monte_carlo_samples_num{50};
  double monte_carlo_3sigma_radius{10.0};

  const point monte_carlo_center{0, 0};

  // Currenlty, we haven't passed these parameters to matlab, so 
    // 1. make them const
    // 2. no variables created for them
  const double model_container_width{200};
  const double model_container_height{200};
  const double model_input_region_height{50};
  const double model_output_region_height{50};

/*------ parameters ------*/

  void dump(std::ostream& os) override;
};

//  Class : GenManager
//-----------------------------------------------------------------------------
class GenManager
{
public:
  GenManager() = delete;

  static Recipe random_recipe(GenSetting* setting);

  static bool 
  launch_matlab(
    GenSetting* setting, 
    const GenDesign &design
  );

  static bool 
  launch_simulation(
    GenSetting* setting, 
    const GenDesign &design
  );

  static std::filesystem::path 
  get_gen_folder(GenSetting* setting, int gid);

  static std::string
  get_design_name(GenSetting* setting, int gid, int design_id);

  static std::filesystem::path
  get_simulation_points_path(GenSetting *setting);

  static std::filesystem::path
  get_recipe_path(const std::filesystem::path & design_folder);

  static std::filesystem::path
  get_velcoity_path(const std::filesystem::path & design_folder);

  static std::filesystem::path
  get_score_path(const std::filesystem::path & design_folder);

  static std::filesystem::path
  get_sim_setting_path(const std::filesystem::path & design_folder);

  static float
  get_simulation_startpoint_z(GenSetting *setting);

  static float
  get_simulation_boundary_z(GenSetting *setting);

  static SimSetting3d
  get_simulation_setting(GenSetting *setting, const GenDesign &design);

  static void
  dump_sim_setting(
    const SimSetting3d & seim_setting, 
    const std::filesystem::path &path
  );

  static float 
  rand_range(float begin, float end);

};

//  Class : GeneticOpt
//-----------------------------------------------------------------------------
class GeneticOpt
{
public:
  GeneticOpt(GenSetting *setting);
  bool run();

  bool restore_session(const std::filesystem::path &path);

  void update_generation(int gid);

  void _evaluate_design(int gid);

  void _gen_next_generation(int gid);

  std::vector<GenDesign*> get_sorted_design(int gid);

  inline int get_generation_num() const;

private:

  GenDesign & _get_design(size_t gid, size_t did);

  void _gen_monte_carlo_samples();

  void _gen_first_generation();

  Recipe _random_merge_two_recipe(const Recipe& r1, const Recipe &r2);

  std::pair<Recipe, Recipe>
  _random_merge_two_recipe_v2(const Recipe& r1, const Recipe &r2);

  void _legalize_recipe(Recipe &recipe);

  void _recipe_mutation(Recipe &recipe, std::string &log);

  void _build_comsol_model(int gid);

  bool _is_generation_folder(const std::filesystem::path &path, int &gid);

  bool _restore_generation_folder(const std::filesystem::path &path, int gid);
  bool _restore_design_folder(
    const std::filesystem::path &design_path, 
    int gid, 
    int design_id
  );
  bool _restore_monte_carlo_samples(const std::filesystem::path &path);

  std::pair<GenDesign*, GenDesign*>
  _randomly_select_two_designs(const std::vector<GenDesign*> &parents);

private:

  GenSetting *_setting;

  int _current_generation_id{0};
  std::vector<std::vector<GenDesign>> _designs;

  std::vector<point3d> _simulating_points;
};

inline int 
GeneticOpt::get_generation_num() const
{
  return _current_generation_id;
}

}

}

#endif
