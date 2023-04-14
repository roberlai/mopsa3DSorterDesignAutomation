#include <mopsa/util/setting.hpp>
#include <mopsa/logger/logger.hpp>

namespace mopsa
{

//  Class : SettingBase
//-----------------------------------------------------------------------------
SettingBase::SettingBase()
  : Reader("%,;[]=\'\"", " \r\t\n")
{

}

bool
SettingBase::read(const std::filesystem::path &path)
{
  LOG(INFO) << "Read setting from " << path << '\n';
  if(!open_file(path)) {
    LOG(ERROR) << "Cannot open " << path << '\n';
    return false;
  }

  bool ok = true;
  std::string token;
  while(!_is_end()) {
    token = read_token();
    if(token == "") break;

    if(auto it = _registered_var.find(token); it != _registered_var.end()) {

      _expect_token(read_token(), "=");
      auto & variable = it->second;
      switch (variable.type) {

        case SettingVarType::Int:
          read_int(* (int*)variable.value_ptr);
          break;

        case SettingVarType::Float:
          read_double(* (double*)variable.value_ptr);
          break;

        case SettingVarType::Boolean:
          read_boolean(* (bool*)variable.value_ptr);
          break;

        case SettingVarType::String:
          read_string(* (std::string*)variable.value_ptr);
          break;

        case SettingVarType::FloatVec:
          _read_vec<double>(variable);
          break;

        case SettingVarType::IntVec:
          _read_vec<int>(variable);
          break;

        case SettingVarType::BooleanVec:
          _read_vec<bool>(variable);
          break;

        default:
          LOG(ERROR) << "Unknow variable type\n";
      } // switch
      _expect_token(read_token(), ";");
    }
    else {
      LOG(WARNING) << "Unknown variable: " << token << '\n';
      while(!_is_end()) {
        token = read_token();
        LOG(WARNING) << "Ignore " << token << '\n';
        if(token == ";") break;
      }
    }
  }

  if(!ok) {
   LOG(ERROR) << "Cannot parse setting: " << path << '\n';
    return false;
  }
  return true;
}

bool 
SettingBase::_is_comment_prefix()
{
  if(_read_cur < _length) {
    return ( _buffer[_read_cur] == '%');
  }
  return false;
};

std::string 
SettingBase::read_comment() 
{
  char c;
  std::string token;
  if(_read_cur < _length) {
    // %
    if(_buffer[_read_cur] == '%') {
      _get_char();
      while(!_is_end()) {
        c = _get_char();
        if(c=='\n') break;
        token.push_back(c);
      }
    } 
  }

  return token;
}

bool 
SettingBase::read_string(std::string &res)
{
  while(!_is_end() and _is_sep_char()) _get_char();
  res = "";
  char c = _cur_char();
  if(c != '\"' and c != '\'') {
    LOG(WARNING) << "Cannot read string at " <<
      _filename << ":" << _line_no << '\n';
    return false;
  }

  bool ok = false;
  char bc = _get_char(); // '\"' or '\''
  while(!_is_end()) {
    c = _get_char();
    if(c == bc) {
      ok = true;
      break;
    }
    else res.push_back(c);
  }

  if(!ok) {
    LOG(WARNING) << "Fail to read string, expect '\"' at " <<
      _filename << ":" << _line_no << '\n';
  }

  return ok;
}


}
