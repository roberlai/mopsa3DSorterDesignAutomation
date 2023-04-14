#include <cstdlib>
#include <iostream>

#include <mopsa/logger/logger.hpp>
#include <mopsa/logger/debug.hpp>
#include <mopsa/chip/flow.hpp>
#include <mopsa/util/util_funcs.h>

namespace mopsa
{

//  Class : FlowBase
//-----------------------------------------------------------------------------

std::string 
FlowBase::read_comment() 
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
FlowBase::_is_comment_prefix()
{
  if(_read_cur < _length) {
    return ( _buffer[_read_cur] == '%');
  }
  return false;
}

//  Class : Flow
//-----------------------------------------------------------------------------

bool
Flow::load_flow(
  const std::filesystem::path &file_x, 
  const std::filesystem::path &file_y
)
{
  open_file(file_x);
  while(!_is_end()) 
  {
    double x, y, v;
    std::string token = read_token();
    if(token.size() == 0) break;

    x = std::stod(token);
    if(!_expect_token(read_token(), ",")) return false;
    y = std::stod(read_token());
    if(!_expect_token(read_token(), ",")) return false;
    v = std::stod(read_token());
    _nodes.emplace_back(point(x, y), v, 0);
  }

  open_file(file_y);
  size_t cnt = 0;
  while(!_is_end()) 
  {
    double x, y, v;
    std::string token = read_token();
    if(token.empty()) break;

    x = std::stod(token);
    if(!_expect_token(read_token(), ",")) return false;
    y = std::stod(read_token());
    if(!_expect_token(read_token(), ",")) return false;
    v = std::stod(read_token());

    ASSERT(_nodes[cnt].coord() == point(x, y));
    _nodes[cnt].vy(v);

    cnt += 1;
  }
  ASSERT(cnt == _nodes.size());

  return true;
}

//  Class : Flow3d
//-----------------------------------------------------------------------------

bool 
Flow3d::load_flow(const std::filesystem::path &file) 
{
  point3d min_pos, max_pos;
  if(!open_file(file)) {
    std::cout << "Cannot open: " << file << std::endl;
  }
  while(!_is_end()) 
  {
    double x, y, z;
    double vx, vy, vz;
    std::string token = read_token();
    if(token.size() == 0) break;

    x = std::stod(token);

    if(!_expect_token(read_token(), ",")) return false;
    y = std::stod(read_token());

    if(!_expect_token(read_token(), ",")) return false;
    z = std::stod(read_token());

    if(!_expect_token(read_token(), ",")) return false;
    vx = std::stod(read_token());

    if(!_expect_token(read_token(), ",")) return false;
    vy = std::stod(read_token());

    if(!_expect_token(read_token(), ",")) return false;
    vz = std::stod(read_token());

    _nodes.emplace_back(point3d(x, y, z), velocity3d(vx, vy, vz));
    const auto& node = _nodes.back();
    if(_nodes.size() == 1) {
      min_pos = max_pos = node.coord();
    }
    else {
      min_pos.x() = std::min(min_pos.x(), node.x());
      min_pos.y() = std::min(min_pos.y(), node.y());
      min_pos.z() = std::min(min_pos.z(), node.z());

      max_pos.x() = std::max(max_pos.x(), node.x());
      max_pos.y() = std::max(max_pos.y(), node.y());
      max_pos.z() = std::max(max_pos.z(), node.z());
    }
  }

  _lower_corner = min_pos;
  _upper_corner = max_pos;

  _length = max_pos.x() - min_pos.x();
  _width  = max_pos.y() - min_pos.y();
  _height = max_pos.z() - min_pos.z();

  Logger::add_record("Flow3d", "length: " + mopsa::to_string(_length), 1);
  Logger::add_record("Flow3d", "width: " + mopsa::to_string(_width), 1);
  Logger::add_record("Flow3d", "height: " + mopsa::to_string(_height), 1);
  Logger::add_record("Flow3d", "lower corner: " + _lower_corner.to_string(), 1);
  Logger::add_record("Flow3d", "upper corner: " + _upper_corner.to_string(), 1);

  return true;
}

}
