# MOPSA3D - simulation

## Simulation Result Example

### With wall collision handling
https://github.com/roberlai/mopsa3DSorterDesignAutomation/assets/56212704/4d414c6c-91aa-4171-9c35-e48747669328

### Without wall collision handling
https://github.com/roberlai/mopsa3DSorterDesignAutomation/assets/56212704/04f1a731-8e69-4d48-8465-af5cdfaceafd

## Compile
1. `mkdir build`
2. `cd build`
3. `cmake ../`
4. `make -j16`
- The binary location: `binary/`.

### System Requirement
- cmake 3.16
- g++ 9.2.1 
- boost 1.65.1
- OpenMp 4.5

### Use docker
  Compile mopsa and simulate chip in docker. Plot result in local machine.
  1. clone project : `git clone https://github.com/roberlai/mopsa.git`
  2. `cd mopsa`
  3. Build docker: `docker build -t mopsa .`
  4. Launch docker and mount mopsa project into docker: `docker run -it -v {your mopsa path}:/mopsa mopsa /bin/bash`
  5. `cd mopsa`
  6. `mkdir build`
  7. `cd build`
  8. `cmake ../`
  9. `make`
  10. Simulate chip: `../scripts/runMopsa.sh ../sim_setting2d/setting_80um.m`
  11. Exit docker: `exit` 
  12. Plot result  `python ../scripts/plot_sim_result.py 80um`

## Setting parameters

- Example 3D simulation setting is under `sim_setting3d`

### Common setting
| Parameter           |  Type  | Default  | Description |
| ------------------- | ------ | -------- | ----------- |
|**chip_name**       | string | "default_name" | Used in output file name |
|**dPs**             | List of float | None | The diameters of particle to simulate |
|**boundary_max_timestep**| int | 100000 | The maximum simulation step|
|**time_resolution**|  float | 1 | Time resolution in simulation|
|**output_folder**| string| None | A path for dumping the simulation result|
|**disable_wall_effect**| boolean | false | Disable wall effect adjustment in simulation |
|**num_threads**| int | 1 | Maximum thread in simulation|
|**show_simulation_progress**| boolean | 0 | Show particle position of each time step (May increase runtime) |

### 3D setting
| Parameter           |  Type  | Default  | Description |
| ------------------- | ------ | -------- | ----------- |
|**mesh_nodes_path** | string |  None | The velocity files of nodes (It should contains velocity at x, y, and z direction.) |
|**init_position_x** | float |  0 | Initial particle position at x direction|
|**init_position_y** | float |  0 | Initial particle position at y direction|
|**init_position_z** | float |  0 | Initial particle position at z direction|
|**boundary_less_z** | float  | None | Simulation boundary. When z position of particle is lesser than **boundary_less_z**, we stop simulating.|

## Run 3D simulation
- set environment variable `MOPSA_ROOT` to the path of the mopsa folder
- `mkdir build`
- `cd build`
- `../scripts/runMopsa3D.sh <setting_path>`
> `../scripts/runMopsa3D.sh ../sim_setting3d/setting_3D.m`

## Plot 3D Simulation Result
- `cd scripts/`
- `python plot_3d_sim_result_vispy.py <sim_output_folder> [num_traj=100]` 
> `python plot_3d_sim_result_vispy.py ../sim_output3d/3D 100`

# MOPSA3D - 3D sorter design automation

## Compile
- Same as the MOPSA3D - simulation

## Setting parameters
- Available parameters are in ./3dsorter_gen_setting/*.m. 

## Before Launch 3D generation
- Before launch 3D sorter generation, we should set the following environment variables
  - `MOPSA_ROOT`: the path to the mopsa folder
  - `COMSOL_ROOT`: the path to comsol folder, something like `comsol55/multiphysics/`
- Make sure `matlab` is runnable directly in the command line

## Launch 3D generation
- `bin/generate3DSorter <setting_path>`
- Run `test_gen_3dsorter_setting.m` first to make sure the flow work.
> `bin/generate3DSorter 3dsorter_gen_setting/test_gen_3dsorter_setting.m`

## Plotting result
- `scripts/plot_ga_score.py`: Plot each generation result
- `scripts/plot_feature_distribution.py`: Plot features distribution

