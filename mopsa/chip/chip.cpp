#include <cstdlib>
#include <mopsa/chip/chip.hpp>
#include <mopsa/logger/logger.hpp>
#include <mopsa/geometry/point.hpp>
#include <mopsa/chip/flow.hpp>
#include <mopsa/headerdef.hpp>

namespace mopsa
{

//  Class : Chip
//-----------------------------------------------------------------------------

Chip::Chip()
{

}

bool 
Chip::load_design(const std::filesystem::path &design_path)
{
  LOG(INFO) << "Load chip design from " << design_path << '\n';
  bool res = _design.load_design(design_path);

  Logger::add_record("chip", "# of obstacles", _design.num_obstacles());
  Logger::add_record("chip", "width:", _design.width());
  Logger::add_record("chip", "height:", _design.height());
  Logger::add_record("chip", "upper left : "
     + _design.upper_left().to_string(), 1);
  Logger::add_record("chip", "lower right : " 
      + _design.lower_right().to_string(), 1);

  return res;
}

bool 
Chip::load_flow(
  const std::filesystem::path &flow_x_path,
  const std::filesystem::path &flow_y_path
)
{
  LOG(INFO) << "Load flow from x: " << flow_x_path 
    << ", y: " << flow_y_path << '\n';

  bool res = _flow.load_flow(flow_x_path, flow_y_path);

  Logger::add_record("chip", "# of flow nodes", _flow.size());

  LOG(INFO) << _flow.size() << " nodes are loaded.\n";
  return res;
}

//  Class : Chip3d
//-----------------------------------------------------------------------------

Chip3d::Chip3d()
{

}

bool 
Chip3d::load_flow(const std::filesystem::path &flow_path)
{
  LOG(INFO) << "Load 3d flow from " << flow_path << '\n';

  bool res = _flow.load_flow(flow_path);

  Logger::add_record("chip", "# of flow nodes", _flow.size());

  LOG(INFO) << _flow.size() << " nodes are loaded.\n";
  return res;
}

}
