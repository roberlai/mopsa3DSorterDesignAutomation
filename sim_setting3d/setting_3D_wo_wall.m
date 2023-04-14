
chip_name = "3D";

mesh_nodes_path = "../data3d/3D.csv";

dPs = [5, 10, 30, 50];

time_resolution = 2;

init_position_x =  0;
init_position_y =  0;
init_position_z = 300;

boundary_max_timestep = 1000000;
boundary_less_z = -5;

disable_wall_effect = true;
output_folder = "../sim_output/3D_wo_wall";

