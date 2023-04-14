#include <mopsa/sim/sim.hpp>

#include <iostream>

int main(int argc, char* argv[])
{
  std::filesystem::path setting_path = argv[1];

  mopsa::Chip3d chip;
  mopsa::SimSetting3d setting3d;

  if(setting3d.read(setting_path)) setting3d.dump(std::cout);
  else mopsa::mopsa_exit(-1);

  chip.load_flow(setting3d.mesh_nodes_path);

  mopsa::Simulate3d(&chip, &setting3d).simulate();
}
