#ifndef MOPSA_UTIL_SETTING_H
#define MOPSA_UTIL_SETTING_H

#include <unordered_map>

#include <mopsa/headerdef.hpp>
#include <mopsa/util/reader.hpp>
#include <mopsa/logger/logger.hpp>

namespace mopsa
{

enum class SettingVarType {
  Int,         /* int */
  Float,       /* double */
  Boolean,     /* bool */
  String,      /* std::string */
  FloatVec,    /* std::vector<double> */
  IntVec,      /* std::vector<int> */
  BooleanVec,  /* std::vector<bool> */
};

class SettingBase : public Reader
{
public:
  SettingBase();

  bool read(const std::filesystem::path &path);
  virtual void dump(std::ostream &os) = 0;

protected:
  struct variable {
    SettingVarType type;
    void *value_ptr;
  };

  template<typename T>
  void _dump_vec(std::ostream &os, std::vector<T>& vec) const;

  template<typename T>
  void _register_variable(
    const std::string &name, 
    SettingVarType type, 
    T &var
  );

private:
  template<typename T>
  bool _read_vec(variable &var);
  /* read string (braced by " or ' (e.g.,  "context" or 'context') ) */
  bool read_string(std::string&) override;

  std::string read_comment() override;

  bool _is_comment_prefix() override;

protected:
  std::unordered_map<std::string, variable> _registered_var;
};

/******************************************************************************
  Inline Function
******************************************************************************/

// Class : SettingBase
// ----------------------------------------------------------------------------
template<typename T>
void 
SettingBase::_register_variable(
  const std::string &name, 
  SettingVarType type, 
  T &var_ptr
)
{
  if(auto it = _registered_var.find(name); it == _registered_var.end()) {

    if constexpr (std::is_same_v<T, int>) {
       ASSERT_MESS(type == SettingVarType::Int, name);
    }
    else if constexpr (std::is_same_v<T, double>) {
      ASSERT_MESS(type == SettingVarType::Float, name);
    }
    else if constexpr (std::is_same_v<T, bool>) {
      ASSERT_MESS(type == SettingVarType::Boolean, name);
    }
    else if constexpr (std::is_same_v<T, std::string>) {
      ASSERT_MESS(type == SettingVarType::String, name);
    }
    else if constexpr (std::is_same_v<T, std::vector<int>>) {
      ASSERT_MESS(type == SettingVarType::IntVec, name);
    }
    else if constexpr (std::is_same_v<T, std::vector<double>>) {
      ASSERT_MESS(type == SettingVarType::FloatVec, name);
    }
    else if constexpr (std::is_same_v<T, std::vector<bool>>) {
      ASSERT_MESS(type == SettingVarType::BooleanVec, name);
    }
    else {
      LOG(DEBUG) << "Invalid variable type: " << name << "\n";
      ASSERT(false);
    }
    _registered_var.insert({name, {type, (void*)&var_ptr}});
  }
  else { 
    LOG(DEBUG) << "Variable: " << name << " has been registered\n";
  }
}

template<typename T>
bool
SettingBase::_read_vec(variable &var) 
{
  bool ok = true;
  _expect_token(read_token(), "[");
  ((std::vector<T>*)(var.value_ptr))->clear();
  while(true) {
    T value;
    if constexpr (std::is_same_v<T, double>) {
      ok = read_double(value);
    } 
    else if constexpr (std::is_same_v<T, int>) {
      ok = read_int(value);
    }
    else if constexpr (std::is_same_v<T, bool>) {
      ok = read_boolean(value);
    }
    else ASSERT(false);

    if(!ok) return false;
    ((std::vector<T>*)(var.value_ptr))->push_back(value);
    
    const auto & token = read_token();
    if(token == "]") break;
    else _expect_token(token, ",");
  }
  return true;
}

template<typename T>
void 
SettingBase::_dump_vec(
  std::ostream &os, 
  std::vector<T>& vec
) const
{
  if(vec.size() == 0) os << "[]";
  else {
    os << "[";
    os << vec.front();
    for(size_t i=1; i<vec.size(); i++)  os << ", " << vec[i];
    os << "]";
  }
}

}

#endif
