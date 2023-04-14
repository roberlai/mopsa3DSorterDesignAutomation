
#include <mopsa/chip/design.hpp>
#include <mopsa/logger/logger.hpp>

#include <iostream>

namespace mopsa
{

/**************************************
 Obstacle
 **************************************/
Obstacle::Obstacle(polygon poly)
  :_polygon(poly)
{

}

/**************************************
 Design
 **************************************/
bool 
Design::load_design(const std::filesystem::path &design_path)
{
  FILE* fp = fopen(design_path.c_str(), "r");

  if(!fp) {
    LOG(ERROR) << "Cannot open " << design_path << " when load design.\n";
    return false;
  }

  int res = 0, num_polygon = 0;
  res = fscanf(fp, "%d", &num_polygon);
  if(res != 1) {
    LOG(ERROR) << "Cannnot read " << design_path << "\n"; 
    return false;
  }

  std::vector<int> a;
  std::vector<double> xs, ys;
  std::vector<polygon> polygons;
  polygons.resize(num_polygon);
  for(int i=0; i<num_polygon; i++) 
  {
    int num_points = 0;
    res = fscanf(fp, "%d", &num_points);
    if(res != 1) {
      LOG(ERROR) << "Cannnot read " << design_path << "\n"; 
      return false;
    }

    polygons[i].reserve(num_points);
    for(int j=0; j<num_points; j++) 
    {
      double x, y;
      res = fscanf(fp, "%lf %lf", &x, &y);
      if(res != 2) {
        LOG(ERROR) << "Cannnot read " << design_path << "\n"; 
        return false;
      }
      polygons[i].push_back({x, y});
      xs.push_back(x);
      ys.push_back(y);
    }

    if(num_points == 1)
    {
      LOG(WARNING) << "Obstacle " << i << " is a point.\n";
    }
    else if(num_points == 2) 
    {
      LOG(WARNING) << "Obstacle " << i << " is a line segment.\n";
    }
    else if(polygons[i].startpoint() != polygons[i].endpoint()) 
    {
      LOG(WARNING) << "Obstacle " << i << " is not a ring. " 
        << polygons[i].startpoint()
        << " <-> "
        << polygons[i].endpoint() << '\n';
    }
  }

  for(auto &poly : polygons) {
    _obstacles.emplace_back(std::move(poly));
  }

  std::sort(xs.begin(), xs.end());
  std::sort(ys.begin(), ys.end());

  _upper_left  = {xs.front(), ys.front()};
  _lower_right = {xs.back(), ys.back()};

  _min_x = xs.front();
  _max_x = xs.back();
  _min_y = ys.front(); 
  _max_y = ys.back();

  _width  = xs.back() - xs.front();
  _height = ys.back() - ys.front();

  LOG(INFO) << "Design starts from " << _upper_left << " to " 
    << _lower_right << '\n';

  return true;
}

}
