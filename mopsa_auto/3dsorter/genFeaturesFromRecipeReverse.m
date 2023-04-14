
function [gen_radius gen_H gen_pos gen_theta gen_phi gen_name] = ...
  genFeaturesFromRecipeReverse()

  global s_feature_container_width
  global s_feature_container_height
  global s_feature_container_radius

  global s_feature_start_x
  global s_feature_start_y
  global s_feature_start_z

  global s_feature_end_x
  global s_feature_end_y
  global s_feature_end_z

  global s_input_region_object
  global s_input_region_width
  global s_input_region_height

  global s_output_region_object
  global s_output_region_width
  global s_output_region_height

  global s_feature_r
  global s_feature_y_init
  global s_feature_y_init_shift
  global s_feature_y_gap

  global s_feature_x_init
  global s_feature_x_gap

  global s_feature_theta
  global s_feature_phi
  global s_feature_fixed_side

  global s_x_factor

  gen_radius = [];
  gen_H = [];
  gen_pos = [];
  gen_theta = [];
  gen_phi = [];
  gen_name = [];

  s_feature_y_gap_p2p = s_feature_y_gap + s_feature_r * 2;
  s_feature_x_gap_p2p = s_feature_x_gap + s_feature_r * 2;
  s_feature_y_init_shift_p2p = s_feature_y_init_shift + s_feature_r * 2;

  s_feature_y_gap_p2p = s_feature_y_gap_p2p / cos(deg2rad(s_feature_phi)) 
  s_feature_x_gap_p2p = s_feature_x_gap_p2p / cos(deg2rad(s_feature_theta - 90)) 

  s_feature_H = sqrt(s_feature_container_width^2 + s_feature_container_width^2) * 1.2;

  s_feature_change_xy = 1
  % randomly choose fixed axis (x, y or z)
  % fixed_side = 1: fixed x
  % fixed_side = 2: fixed y
  % fixed_side = 3: fixed z

  c = mod(s_feature_fixed_side, 3);
  c1 = c + 1 
  c2 = mod(c+1, 3) + 1
  c3 = mod(c+2, 3) + 1 %fixed axis

  space = s_feature_r * s_x_factor * 1.5;
  i = 0;
  x0 = s_feature_end_x - s_feature_r
  for x=x0+space:-s_feature_x_gap_p2p:s_feature_start_x-space
    j = 0;
    y0 = s_feature_start_y + ... 
      mod(s_feature_y_init + s_feature_y_init_shift_p2p * i, s_feature_y_gap_p2p) - space;
    for y=y0:s_feature_y_gap_p2p:s_feature_end_y+space

      %if s_feature_fixed_side > 2
        %z = s_feature_start_z;
      %else
        %z = s_feature_end_z;
      %end

      z = (s_feature_end_z + s_feature_start_z)/2

      pos = [0 0 0];
      if s_feature_change_xy == 1
        pos(c1) = y;
        pos(c2) = x;
        pos(c3) = z;
      else
        pos(c1) = x;
        pos(c2) = y;
        pos(c3) = z;
      end

      name = "feature" + string(i) + "_" + string(j);
      gen_pos = [gen_pos; pos];
      gen_name = [gen_name name];
      gen_radius = [gen_radius s_feature_r];
      gen_H = [gen_H s_feature_H];
      gen_theta = [gen_theta s_feature_theta];
      gen_phi = [gen_phi s_feature_phi];

      name = "feature" + string(i) + "_" + string(j) + "_r";
      gen_pos = [gen_pos; pos];
      gen_name = [gen_name name];
      gen_radius = [gen_radius s_feature_r];
      gen_H = [gen_H s_feature_H];
      gen_theta = [gen_theta s_feature_theta + 180];
      gen_phi = [gen_phi s_feature_phi];

      j = j + 1;
    end
    i = i + 1;
  end

end



