#include <mopsa/gen/gen_3dsorter.hpp>
#include <mopsa/gen/gen_3design.hpp>
#include <mopsa/logger/logger.hpp>

#include <random>

using namespace mopsa::Gen3DSorter;

namespace mopsa::Gen3DSorter {

std::ostream& 
operator << (std::ostream &out, const GenDesignStatus& status)
{
  switch (status) {
    case GenDesignStatus::UNINIT:
      out << "UNINIT";
      break;

    case GenDesignStatus::BUILT_RECIPE:
      out << "BUILT_RECIPE";
      break;

    case GenDesignStatus::BUILT_COMSOL:
      out << "BUILT_COMSOL";
      break;

    case GenDesignStatus::BUILT_COMSOL_FAILED:
      out << "BUILT_COMSOL_FAILED";
      break;

    case GenDesignStatus::SIMULATED:
      out << "SIMULATED";
      break;

    default:
      out << "UNINIT";
  }

  return out;
}

}

//  Class : Recipe
//-----------------------------------------------------------------------------
void 
Recipe::dump_recipe(
  const std::filesystem::path &output_path
) const
{
  FILE *fp = fopen(output_path.c_str(), "w");

  fprintf(fp, "radius, theta, phi,"
    "y_init, y_init_shift, y_gap, x_init, x_gap, fixed_side\n"
  );
  fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d\n",
    radius, theta, phi,
    y_init, y_init_shift, y_gap,
    x_init, x_gap, 
    fixed_side
  );
  fclose(fp);
}

bool
Recipe::read_recipe(
  const std::filesystem::path &path
)
{
  FILE *fp = fopen(path.c_str(), "r");
  LOG(INFO) << "Load recipe from " << path << '\n';

  char title[200];
  ASSERT(fgets(title, 200, fp) != nullptr);

  ASSERT(fscanf(fp, "%f %f %f %f %f %f %f %f %d\n",
    &radius, &theta, &phi,
    &y_init, &y_init_shift, &y_gap,
    &x_init, &x_gap, 
    &fixed_side
  ) == 9);

  fclose(fp);
  return true;
}

//  Class : GenDesign
//-----------------------------------------------------------------------------

GenDesign::GenDesign()
{

}

void 
GenDesign::init(
  int gen_id,
  int design_id,
  const std::string &design_name,
  const std::filesystem::path &design_folder,
  const std::filesystem::path &recipe_path,
  const Recipe& recipe
)
{
  _gen_id = gen_id;
  _design_id = design_id;
  _design_name = design_name;
  _design_folder = design_folder;
  _recipe = recipe;
  _recipe_path = recipe_path;
  _status = GenDesignStatus::BUILT_RECIPE;

  std::filesystem::create_directory(design_folder);
  recipe.dump_recipe(recipe_path);
}

bool
GenDesign::restore_session(
  int gen_id,
  int design_id,
  const std::string &design_name,
  const std::filesystem::path &design_folder
)
{
  _gen_id = gen_id;
  _design_id = design_id;
  _design_name = design_name;
  _design_folder = design_folder;

  const auto & recipe_path = GenManager::get_recipe_path(design_folder);
  if(std::filesystem::exists(recipe_path) 
    && _recipe.read_recipe(recipe_path)) 
  {
    _recipe_path = recipe_path;
    _status = GenDesignStatus::BUILT_RECIPE;
  }

  if(_status == GenDesignStatus::BUILT_RECIPE) {
    update_vel_path();
  }

  if(_status == GenDesignStatus::BUILT_COMSOL) {
    update_score();
  }

  return true;
}

bool
GenDesign::update_vel_path()
{
  const auto & vel_path = GenManager::get_velcoity_path(_design_folder);
  if(_status == GenDesignStatus::BUILT_RECIPE) 
  {
    if(std::filesystem::exists(vel_path)) {
      _velocity_nodes_path = vel_path;
      _status = GenDesignStatus::BUILT_COMSOL;
    }
    else if(std::filesystem::exists(_design_folder / "comsol_fail")) {
      _status = GenDesignStatus::BUILT_COMSOL_FAILED;
    }
  }

  return true;
}

bool
GenDesign::update_score()
{
  const auto & score_path = GenManager::get_score_path(_design_folder);
  if(_status == GenDesignStatus::BUILT_COMSOL
    && std::filesystem::exists(score_path)
  ) 
  {
    _score_path = score_path;
    _status = GenDesignStatus::SIMULATED;

    FILE * fp = fopen(score_path.c_str(), "r");
    ASSERT(fscanf(fp, "%f\n", &_score) == 1);
    fclose(fp);

    if(_score >= 10000) _score = 0;

    return true;
  }

  return false;
}

void
GenDesign::copy_result_from(const GenDesign& design)
{
  _score = design._score;

  _velocity_nodes_path = GenManager::get_velcoity_path(_design_folder);
  std::filesystem::copy(
    design.get_velcoity_nodes_path(), 
    _velocity_nodes_path,
    std::filesystem::copy_options::overwrite_existing
  );

  _score_path = GenManager::get_score_path(_design_folder);

  std::filesystem::copy(
    std::filesystem::absolute(design.get_design_folder() / "sim_output"),
    std::filesystem::absolute(_design_folder / "sim_output"),
    std::filesystem::copy_options::overwrite_existing  |
    std::filesystem::copy_options::recursive
  );

  _status = design._status;

  const auto & evolution_log = _design_folder / "evolution.log";
  FILE *fp = fopen(evolution_log.c_str(), "w");
  fprintf(fp, "copy from: %s\n", design.get_design_name().c_str());
  fclose(fp);
}

