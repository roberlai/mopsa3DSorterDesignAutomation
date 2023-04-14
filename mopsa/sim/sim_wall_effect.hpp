#include <omp.h>
#include <mopsa/sim/sim.hpp>
#include <vector>

namespace mopsa
{

enum class WallEffectStatus
{
  uninit = 0,
  has_collision,
  no_collision,
  time_resolution_too_large 
};

static std::vector<std::string> WallEffectStatusString =
{
  "uninit",
  "has collision",
  "no collision",
  "time resolution is too large"
};

//  Class: GeneralWallEffect
//  Abstract: Perform general wall effect in a independent class
//-----------------------------------------------------------------------------
template<class P, class G, class F>
class GeneralWallEffect
{

public:

  using point            = typename F::point_type;
  using flow_grid        = G;
  using flow_block_group = typename flow_grid::flow_block_group_type;

  GeneralWallEffect(
    const P & particle, 
    const G & flow_grid, 
    const F & flow,
    flow_block_group & group,
    double collision_distance_threshold,
    int used_thread,
    bool debug=false
  );

  bool run();

  // query data
  inline const auto & group() const;
  inline const auto & num_on_surface_nodes() const;
  inline const auto & num_in_particle_nodes() const;
  inline const auto & distance_adjustment() const;
  inline const auto & adjsut_direction() const;
  inline const auto & new_coord() const;
  inline bool has_collision() const;
  inline WallEffectStatus status() const;
  inline std::string get_status_string() const;

private:

  inline void _init();

private:

  const P & _particle;
  const G & _flow_grid;
  const F & _flow;
  flow_block_group & _group;

  double _collision_distance_threshold;
  //double _distance_to_adjustment();

  int _used_thread{1};

  int _on_surface_sum;
  int _in_particle_sum;

  bool _has_collision;
  double _adj_dist;
  point _adj_dir;
  point _new_coord;

  WallEffectStatus _status;

  bool _debug{0};
};

/******************************************************************************
  Inline Function
******************************************************************************/

//  Class : GeneralWallEffect
//-----------------------------------------------------------------------------

template<class P, class G, class F>
GeneralWallEffect<P, G, F>::GeneralWallEffect(
  const P & particle, 
  const G & flow_grid,
  const F & flow,
  flow_block_group & group,
  double collision_distance_threshold,
  int used_thread,
  bool debug
)
  : _particle(particle)
  , _flow_grid(flow_grid)
  , _flow(flow)
  , _group(group)
  , _collision_distance_threshold(collision_distance_threshold)
  , _used_thread(used_thread)
  , _debug(debug)
{

  if(debug) {
    printf("In %s, collision_distance_threshold: %lf\n",
      __FUNCTION__,
      collision_distance_threshold
    );
  }
}

// Abstract: initial calculation data before calculate wall effect
template<class P, class G, class F>
inline void
GeneralWallEffect<P, G, F>::_init()
{
  _on_surface_sum = 0;
  _in_particle_sum = 0;

  _has_collision = false;

  _adj_dist = 0;
  _adj_dir.set_to_init();
  _new_coord.set_to_init();

  _status = WallEffectStatus::uninit;
}

template<class P, class G, class F>
bool
GeneralWallEffect<P, G, F>::run()
{
  _init();

  _flow_grid.get_adjacent_blocks(&_particle, _group);

  int nodes_sum = _group.size(FlowBlockNodeType::obstacle_node);

  std::vector<std::vector<int>> 
    on_surface_id_t(_used_thread, std::vector<int>());
  std::vector<std::vector<int>> 
    in_particle_id_t(_used_thread, std::vector<int>());
  std::vector<point> on_surface_vec_t(_used_thread);
  std::vector<double> dist_to_surf_t(_used_thread, 0);

  if(_debug) {
    printf("# of obstacle nodes: %d\n", nodes_sum);
  }

  #pragma omp parallel for num_threads(_used_thread)
  for(int i=0; i<nodes_sum; i++) {
    int tid     = _used_thread > 1? omp_get_thread_num():0;
    int node_id = _group.get(FlowBlockNodeType::obstacle_node, i);

    const auto & node = _flow.nodes()[node_id];

    bool on_surface = false, in_particle = false;
    double distance = _particle.distance(node.coord());
    double diff = distance - _particle.radius();

    if(std::abs(diff) < _collision_distance_threshold) {
      on_surface = true;
    } 
    else if (diff <= 0) {
      in_particle = true;
    }

    ASSERT( ! (on_surface && in_particle) );

    if(on_surface) {
      on_surface_id_t[tid].push_back(node_id);
      on_surface_vec_t[tid] += _particle.coord() - node.coord();
    }

    if(in_particle) {
      in_particle_id_t[tid].push_back(node_id);
    }

    if(in_particle) ASSERT(distance <= _particle.radius());
    if(distance <= _particle.radius()) {
      double dist_to_surf = _particle.radius() - distance;
      dist_to_surf_t[tid] = std::max(dist_to_surf_t[tid], dist_to_surf);

      if(_debug) {
        printf("\tNode %d, on_surface: %d, in_particle: %d, dist=%.9lf\n",
          node_id, on_surface, in_particle, distance
        );
        printf("\t_dist_to_surf_t[%d] = %.9lf\n",
          tid, dist_to_surf
        );
      }
    }
  }

  _on_surface_sum = on_surface_id_t[0].size();
  _in_particle_sum = in_particle_id_t[0].size();
  _adj_dist = dist_to_surf_t[0];
  for(int i=1; i<_used_thread; i++) {
    _on_surface_sum += on_surface_id_t[i].size();
    on_surface_vec_t[0] += on_surface_vec_t[i];

    _in_particle_sum += in_particle_id_t[i].size();

    _adj_dist = std::max(_adj_dist, dist_to_surf_t[i]);
  }

  if(_debug) {
    printf("adj_dist: %.9lf\n", _adj_dist);
  }
  bool ok = true;
  _has_collision = false;
  _new_coord = _particle.coord();
  if(_on_surface_sum == 0 and _in_particle_sum == 0) {
    // no collision
    _new_coord = _particle.coord();
    _status = WallEffectStatus::no_collision;
  }
  else if(_on_surface_sum == 0 and _in_particle_sum > 0) {
    //ASSERT(false && "time resolution is too large");
    // TODO: time resolution is too large.
    _new_coord = _particle.coord();
    _status = WallEffectStatus::time_resolution_too_large;
    ok = false;
  }
  else { // on_surface_sum > 0
    ASSERT(_on_surface_sum > 0);

    if(_adj_dist > 0) {
      _has_collision = true;
      _adj_dir = (on_surface_vec_t[0] / _on_surface_sum);
      _adj_dir /= _adj_dir.len();

      auto adj_vec = _adj_dir * _adj_dist;
      if(_debug) {
        std::cout << "Adj vec: " << adj_vec << '\n';
      }
      _new_coord = _particle.coord() + adj_vec;

      _status = WallEffectStatus::has_collision;
    }
    else _status = WallEffectStatus::no_collision;
  }

  return ok;
}

template<class P, class G, class F>
inline const auto & 
GeneralWallEffect<P, G, F>::group() const
{
  return _group;
}

template<class P, class G, class F>
inline const auto &
GeneralWallEffect<P, G, F>::num_on_surface_nodes() const
{
  return _on_surface_sum;
}

template<class P, class G, class F>
inline const auto &
GeneralWallEffect<P, G, F>::num_in_particle_nodes() const
{
  return _in_particle_sum;
}

template<class P, class G, class F>
inline const auto & 
GeneralWallEffect<P, G, F>::distance_adjustment() const
{
  return _adj_dist;
}

template<class P, class G, class F>
inline const auto & 
GeneralWallEffect<P, G, F>::adjsut_direction() const
{
  return _adj_dir;
}

template<class P, class G, class F>
inline const auto & 
GeneralWallEffect<P, G, F>::new_coord() const
{
  return _new_coord;
}

template<class P, class G, class F>
inline bool 
GeneralWallEffect<P, G, F>::has_collision() const
{
  return _has_collision;
}
 
template<class P, class G, class F>
inline WallEffectStatus 
GeneralWallEffect<P, G, F>::status() const
{
  return _status;
}

template<class P, class G, class F>
inline std::string 
GeneralWallEffect<P, G, F>::get_status_string() const
{
  return WallEffectStatusString[static_cast<int>(_status)];
}

};
