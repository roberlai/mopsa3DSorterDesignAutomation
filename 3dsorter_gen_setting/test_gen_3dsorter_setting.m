
output_folder = "test_gen_3dsorter_5_15um";

target_particle_dims = [5, 15];

ga_num_elite = 100;

%ga_initial_population_size = 500;
%ga_population_size = 1000;
%ga_num_iterations = 10;

ga_initial_population_size = 2;
ga_population_size = 2;
ga_num_iterations = 1;

recipe_radius_range_min_in_um = 50;
recipe_radius_range_max_in_um = 60;
recipe_theta_range_min_in_degree = 30;
recipe_theta_range_max_in_degree = 150;
recipe_phi_range_min_in_degree = 0;
recipe_phi_range_max_in_degree = 40;

sim_samples_count = 50;
sim_3sigma_radius = 10;
