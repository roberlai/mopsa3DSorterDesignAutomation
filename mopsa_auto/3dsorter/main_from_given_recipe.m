
link_comsol

import com.comsol.model.*
import com.comsol.model.util.*

num_total = 1
for iter = 1: num_total

  if ~isfile(s_recipe_path)
    fprintf(1, "Cannot open %s\n", s_recipe_path)
    continue;
  end

  s_feature_container_width = 200 * s_x_factor;
  s_feature_container_height = 200 * s_x_factor;
  s_feature_container_radius = s_feature_container_width / 2;

  s_feature_start_x = -s_feature_container_radius;
  s_feature_start_y = -s_feature_container_radius;
  s_feature_start_z = -s_feature_container_radius;

  s_feature_end_x = s_feature_container_radius;
  s_feature_end_y = s_feature_container_radius;
  s_feature_end_z = s_feature_container_radius;

  s_input_region_object = "input_region";
  s_input_region_width = s_feature_container_width;
  s_input_region_height = 50 * s_x_factor;

  s_output_region_object = "output_region";
  s_output_region_width = s_feature_container_width;
  s_output_region_height = 50 * s_x_factor;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mkdir(s_model_path);

  fprintf(1, "Design is created at %s\n", s_model_path)

  [s_feature_r s_feature_y_init s_feature_y_init_shift s_feature_y_gap ...
    s_feature_x_init s_feature_x_gap ...
    s_feature_theta s_feature_phi s_feature_fixed_side] ...
  = readRecipe(s_recipe_path);

  [gen_radius gen_H gen_pos gen_theta gen_phi gen_name] = ...
    genFeaturesFromRecipeReverse();

  s_num_features = length(gen_name);

  s_output_features_recipe_path = s_model_path + "/feature_recipe_from_matlab.txt"

  out = fopen(s_output_features_recipe_path, "w");

  fprintf(out, "r theta phi y_init, y_init_shift, y_gap x_init x_gap fixed_side\n")
  fprintf(out, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %d\n", ...
    s_feature_r, s_feature_theta, s_feature_phi, ...
    s_feature_y_init, s_feature_y_init_shift, s_feature_y_gap, ...
    s_feature_x_init, s_feature_x_gap, ...
    s_feature_fixed_side ...
  );

  fprintf(out, "%d\n", s_num_features);
  for i=1:s_num_features
    fprintf(out, "%s\n", gen_name(i));
    fprintf(out, "%.9f\n", gen_radius(i));
    fprintf(out, "%.9f\n", gen_H(i));
    fprintf(out, "%.9f %.9f %.9f\n", gen_pos(i, 1), gen_pos(i, 2), gen_pos(i, 3));
    fprintf(out, "%.9f\n", gen_theta(i));
    fprintf(out, "%.9f\n", gen_phi(i));
  end

  fclose(out);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %Create model
  ModelUtil.remove(s_model_name);

  model = ModelUtil.create(s_model_name);
  model.modelPath(s_model_path);
  model.component.create('comp1', true);
  model.component('comp1').geom.create('geom1', 3);

  fprintf(1, "Create features\n")
  createFeatures

  createBody

  model.component('comp1').geom('geom1').create('sel1', 'ExplicitSelection');
  model.component('comp1').geom('geom1').feature('sel1').selection('selection').init(2);
  model.component('comp1').geom('geom1').feature('sel1').label('input');
  model.component('comp1').geom('geom1').feature('sel1').selection('selection').set(s_input_region_object, 4);
  model.component('comp1').geom('geom1').run('sel1');

  model.component('comp1').geom('geom1').create('sel2', 'ExplicitSelection');
  model.component('comp1').geom('geom1').feature('sel2').label('output');
  model.component('comp1').geom('geom1').feature('sel2').selection('selection').init(2);
  model.component('comp1').geom('geom1').feature('sel2').selection('selection').set(s_output_region_object, 1);
  model.component('comp1').geom('geom1').run('sel2');

  model.component('comp1').geom('geom1').create('dif1', 'Difference')
  model.component('comp1').geom('geom1').feature('dif1').selection('input').set("body");
  model.component('comp1').geom('geom1').feature('dif1').selection('input2').set(all_features)
  model.component('comp1').geom('geom1').run

  addWaterAsMaterial
  %
  model.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
  model.component('comp1').physics('spf').create('inl1', 'InletBoundary', 2);
  model.component('comp1').physics('spf').feature('inl1').selection.named('geom1_sel1');
  model.component('comp1').physics('spf').create('out1', 'OutletBoundary', 2);
  model.component('comp1').physics('spf').feature('out1').selection.named('geom1_sel2');
  model.component('comp1').physics('spf').feature('inl1').set('U0in', 0.001);
  model.component('comp1').physics('spf').prop('PhysicalModelProperty').set('StokesFlowProp', true);

  %
  s_study_succeeded = 0

  try

    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').automatic(true);
    model.component('comp1').mesh('mesh1').autoMeshSize(s_mesh_size);

    if s_skip_mesh == 0
      model.component('comp1').mesh('mesh1').run;
    else
      fprintf(1, "Skip mesh")
    end

    s_mesh_nodes_filename = "vel_size_" + num2str(s_mesh_size) + ".csv"
    studyAndGetResult

    if s_skip_mesh == 0 && s_skip_study == 0
      model.sol('sol1').runAll;
      model.result.export('data1').run;
      s_study_succeeded = 1
    else
      fprintf(1, "Skip study")
    end

  catch e
    fprintf(1,'%s\n',e.identifier);
    fprintf(1,'%s\n',e.message);
    s_study_succeeded = 0
  end

  if s_study_succeeded == 1 && s_skip_study == 0 && s_skip_mesh == 0
    s_output_image = s_model_path + "/figure.png"
    createImage

    out = fopen(s_model_path + "/comsol_done", "w")
    fprintf(out, "done")
    fclose(out)
  else
    out = fopen(s_model_path + "/comsol_fail", "w")
    fprintf(out, "fail")
    fclose(out)
  end

  fprintf(1, "s_save_model: %d ", s_save_model)
  if s_save_model == 1
    output_path = s_model_path + "/" + s_model_name
    fprintf(1, "save model to " + output_path)
    mphsave(model, output_path)
  end
    
end
