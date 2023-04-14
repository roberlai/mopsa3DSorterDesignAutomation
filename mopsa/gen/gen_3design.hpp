#ifndef MOPSA_GEN_3DDESIGN_H
#define MOPSA_GEN_3DDESIGN_H

#include <filesystem>

namespace mopsa {

namespace Gen3DSorter {

class GeneticOpt;

//  Class : Recipe
//-----------------------------------------------------------------------------
struct Recipe
{
  float radius;
  float theta;
  float phi;

  float y_init;
  float y_init_shift;
  float y_gap;

  float x_init;
  float x_gap;

  int fixed_side;

  void dump_recipe(const std::filesystem::path &output_path) const;

  bool read_recipe(const std::filesystem::path &path);
};

//  Enum Class : GenDesignStatus
//-----------------------------------------------------------------------------
enum class GenDesignStatus
{
  UNINIT,
  BUILT_RECIPE,
  BUILT_COMSOL,
  BUILT_COMSOL_FAILED,
  SIMULATED,
  SIMULATED_FAILED,
};

extern
std::ostream& operator << (std::ostream &out, const GenDesignStatus& status);

//  Class : GenDesign
//-----------------------------------------------------------------------------
class GenDesign
{
public:

  GenDesign();

  void init(
    int gen_id,
    int deisign_id,
    const std::string &design_name,
    const std::filesystem::path &design_folder,
    const std::filesystem::path &recpie_path,
    const Recipe& recipe
  );

  bool restore_session(
    int gen_id,
    int design_id,
    const std::string &design,
    const std::filesystem::path &design_folder
  );

  bool update_vel_path();

  bool update_score();

  void copy_result_from(const GenDesign& design);

  inline const std::string& get_design_name() const;
  inline const std::filesystem::path &get_design_folder() const;
  inline const std::filesystem::path &get_recipe_path() const;
  inline const std::filesystem::path &get_velcoity_nodes_path() const;
  inline const Recipe & get_recipe() const;
  // inline Recipe & get_recipe();
  inline bool is_uninit() const;
  inline GenDesignStatus get_status() const;
  inline float get_score() const;

private:
  GenDesignStatus _status{GenDesignStatus::UNINIT};
  int _gen_id{-1};
  int _design_id{-1};
  std::string _design_name;
  std::filesystem::path _design_folder;
  Recipe _recipe;

  float _score{0};
  std::filesystem::path _recipe_path{""};
  std::filesystem::path _velocity_nodes_path{""};
  std::filesystem::path _score_path{""};
};

inline const std::string& 
GenDesign::get_design_name() const
{
  return _design_name;
}

inline const std::filesystem::path &
GenDesign::get_design_folder() const
{
  return _design_folder;
}

inline const Recipe & 
GenDesign::get_recipe() const
{
  return _recipe;
}

//inline Recipe & 
//GenDesign::get_recipe()
//{
  //return _recipe;
//}

inline const std::filesystem::path&
GenDesign::get_recipe_path() const
{
  return _recipe_path;
}

inline const std::filesystem::path&
GenDesign::get_velcoity_nodes_path() const
{
  return _velocity_nodes_path;
}

inline bool 
GenDesign::is_uninit() const
{
  return _status == GenDesignStatus::UNINIT;
}

inline GenDesignStatus 
GenDesign::get_status() const
{
  return _status;
}

inline float 
GenDesign::get_score() const
{
  return _score;
}

}
  
}
#endif
