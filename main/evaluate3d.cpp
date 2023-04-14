#include <mopsa/sim/sim.hpp>

#include <iostream>

using namespace mopsa;
int main(int argc, char* argv[])
{
  std::filesystem::path setting_path = argv[1];

  mopsa::Chip3d chip;
  mopsa::SimSetting3d setting3d;

  if(setting3d.read(setting_path)) setting3d.dump(std::cout);
  else mopsa::mopsa_exit(-1);

  chip.load_flow(setting3d.mesh_nodes_path);

  mopsa::Simulate3d sim(&chip, &setting3d);
  bool sim_ok = sim.simulate();
  ASSERT(sim_ok);

  double INF = std::numeric_limits<double>().max();
  double score = INF;

  LOG(INFO) << "\n====== Start evaluating ======\n";

  std::map<float, point3d> average_positions;

  const auto & dps = setting3d.dPs;
  const auto & init_pos = sim.get_initial_position();
  for(size_t i = 0; i < dps.size(); i++) {

    float dp = dps[i];

    int fail_cnt = 0;

    point3d ave_pos{0,0,0};
    int suc_cnt = 0;
    for(size_t j = 0; j < init_pos.size(); j++) {
      auto [ok, sim_result] = sim.get_sim_result(j, dp);

      if(sim_result.status == SimStatus::ArrivalExpectedBundary) {
        ave_pos += sim_result.final_position;
        suc_cnt += 1;
      }
      else {
        fail_cnt += 1;
        LOG(INFO) << "Failed: int_pos_id: " << j << " " << init_pos[j]
          << " diam = "<< dp << ", final position is at " 
          << sim_result.final_position 
        << '\n';
      }
    }

    ave_pos /= suc_cnt;
    LOG(INFO) << "Diam " << dp << ", average position: " 
      << ave_pos << " (" << suc_cnt << "/" << init_pos.size() << ")\n";

    if(fail_cnt <= 0.2 * init_pos.size()) {
      average_positions[dp] = ave_pos;
    }
    else {
      average_positions[dp] = {INF, INF, INF};
      LOG(INFO) << "Ignore diam: " << dp << " in evalulation since failed"
        " count: " << fail_cnt << " is more than 20%\n";
    }
  }

  for(size_t i=0; i<init_pos.size(); i++) {

    double min_dist = std::numeric_limits<double>::max();
    for(size_t j=0; j<dps.size(); j++) {
      for(size_t k=j + 1; k<dps.size(); k++) {
        auto [okj, sim_result_j] = sim.get_sim_result(i, dps[j]);
        auto [okk, sim_result_k] = sim.get_sim_result(i, dps[k]);

        if(!okj || !okk) continue;
        
        double dist = boost::geometry::distance(
          sim_result_j.final_position, 
          sim_result_k.final_position
        );

        min_dist =  std::min(dist, min_dist);
      }
    }
    LOG(INFO) << "Pos " << i << " " << init_pos[i] 
      << " min_dist: " << min_dist << '\n';
  }
  LOG(INFO) << '\n';

  for(size_t i=0; i<dps.size(); i++) {
    const auto & pos1 = average_positions[dps[i]];
    if(pos1 == point3d(INF, INF, INF)) {
      score = INF;
      break;
    }
    for(size_t j=0; j<dps.size(); j++) {
      if(i==j) continue;
      const auto & pos2 = average_positions[dps[j]];
      double dist = boost::geometry::distance(pos1, pos2);
      std::cout << dist << std::endl;
      if(dist < score) {
        score = dist;
        LOG(INFO) << "Neaseet pair: " << dps[i] << " " << dps[j] 
          << " dist = " << dist << '\n';
      }
    }
  }
  if(score == INF) score = 0;
  LOG(INFO) << "Final score: " << score << '\n';

  std::filesystem::path output_score_path = setting3d.output_folder;
  output_score_path /= "score.txt";

  LOG(INFO) << "Output score to " << output_score_path << '\n';
  FILE* fp = fopen(output_score_path.c_str(), "w");
  fprintf(fp, "%.3lf\n", score);

  for(const auto &dp: dps) {
    mopsa::point3d pos = average_positions[dp];
    fprintf(fp, "%.2f: %.3lf %.3lf %.3lf\n",
      dp, pos.x(), pos.y(), pos.z()
    );
  }
  fclose(fp);

  return 0;
}
