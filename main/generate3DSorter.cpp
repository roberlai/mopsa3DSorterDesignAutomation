#include <iostream>
#include <cstdlib>
#include <mopsa/sim/sim.hpp>
#include <mopsa/gen/gen_3dsorter.hpp>
#include <filesystem>

using namespace mopsa;

int main(int argc, char* argv[])
{
  std::filesystem::path setting_path = argv[1];

  char * root = getenv("MOPSA_ROOT");
  if(!root) {
    LOG(ERROR) << "Cannot find environment variable: MOPSA_ROOT. "
      "Please set environment variable `MOPSA_ROOT` to the MOPSA folder\n";
    mopsa::mopsa_exit(-1);
  }
  LOG(INFO) << "MOPSA_ROOT: " << root << '\n';

  char * comsol_root = getenv("COMSOL_ROOT");
  if(!comsol_root) {
    LOG(ERROR) << "Cannot find environment variable: COMSOL_ROOT. "
      "Please set environment variable `COMSOL_ROOT` to the COMSOL folder\n";
    mopsa::mopsa_exit(-1);
  }
  LOG(INFO) << "COMSOL_ROOT: " << comsol_root << '\n';

  LOG(INFO) << "Search matlab ...\n";
  if(system("command -v matlab")!=0) {
    LOG(ERROR) << "Cannot find matlab\n";
    mopsa::mopsa_exit(-1);
  }

  LOG(INFO) << "Search comsol ...\n";
  if(system("command -v comsol")!=0) {
    LOG(ERROR) << "Cannot find comsol\n";
    mopsa::mopsa_exit(-1);
  }

  std::filesystem::path mopsa_root(root);
  mopsa::Gen3DSorter::GenSetting setting;
  if(setting.read(setting_path)) {
    setting.dump(std::cout);
  }
  else mopsa::mopsa_exit(-1);

  mopsa::Gen3DSorter::GeneticOpt opt(&setting);
  
  if(std::filesystem::exists(setting.gen_path)) {
    opt.restore_session(setting.gen_path);
    int a;
    printf("press any key to continue.");
    std::cin >> a;
  }

  return opt.run();
}
