
name = 'body'
model.component('comp1').geom('geom1').create(name, 'Block');

model.component('comp1').geom('geom1').feature(name).set('size', [ ...
  s_feature_container_width s_feature_container_width s_feature_container_height]);

model.component('comp1').geom('geom1').feature(name).set('base', 'center');
model.component('comp1').geom('geom1').feature(name).set('pos', [0 0 0]);
model.component('comp1').geom('geom1').run(name);

name = s_input_region_object
model.component('comp1').geom('geom1').create(name, 'Block')
model.component('comp1').geom('geom1').feature(name).label(name);
model.component('comp1').geom('geom1').feature(name).set('size', [...
  s_input_region_width s_input_region_width s_input_region_height]);
model.component('comp1').geom('geom1').feature(name).set('base', 'center');
pos = [0 0 s_feature_container_height/2 + s_input_region_height/2]
model.component('comp1').geom('geom1').feature(name).set('pos', pos);
model.component('comp1').geom('geom1').run(name);

name = s_output_region_object
model.component('comp1').geom('geom1').create(name, 'Block')
model.component('comp1').geom('geom1').feature(name).label(name);
model.component('comp1').geom('geom1').feature(name).set('size', [...
  s_input_region_width s_input_region_width s_input_region_height]);
model.component('comp1').geom('geom1').feature(name).set('base', 'center');
pos = [0 0 -(s_feature_container_height/2 + s_input_region_height/2)]
model.component('comp1').geom('geom1').feature(name).set('pos', pos);
model.component('comp1').geom('geom1').run(name);

% change length unit from m to um
model.geom('geom1').lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
