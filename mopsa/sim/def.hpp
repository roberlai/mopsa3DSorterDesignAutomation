#ifndef MOPSA_SIM_DEF_H
#define MOPSA_SIM_DEF_H

#include <mopsa/headerdef.hpp>
#include <mopsa/sim/sim_flowgrid.hpp>

namespace mopsa
{

template<class T>
class FlowBlock;

template<class T>
class FlowBlockGroup;

using flow_block = FlowBlock<point>;
using flow_block_group = FlowBlockGroup<point>;

using flow_block3d = FlowBlock<point3d>;
using flow_block_group3d = FlowBlockGroup<point3d>;

}

#endif
