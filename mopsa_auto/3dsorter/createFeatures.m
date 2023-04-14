
for i = 1:s_num_features

  name = gen_name(i)
  model.component('comp1').geom('geom1').create(name, 'Cylinder');
  model.component('comp1').geom('geom1').feature(name).set('axistype', 'spherical');

  model.component('comp1').geom('geom1').feature(name).set('r', gen_radius(i));
  model.component('comp1').geom('geom1').feature(name).set('h', gen_H(i));
  model.component('comp1').geom('geom1').feature(name).set('pos', gen_pos(i, :));
  model.component('comp1').geom('geom1').feature(name).set('ax2', [gen_theta(i) gen_phi(i)]);
end

model.component('comp1').geom('geom1').run(name);

all_features = gen_name
