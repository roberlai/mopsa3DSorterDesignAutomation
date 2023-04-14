
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('data', 'dset1');
model.result.export.create('data1', 'Data');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
%model.sol('sol1').runAll;

model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result.export('data1').set('expr', {'comp1.u' 'comp1.v' 'comp1.w'});
model.result.export('data1').set('unit', {'m/s' 'm/s' 'm/s'});
model.result.export('data1').set('descr', {[native2unicode(hex2dec({'90' '1f'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode')  native2unicode(hex2dec({'58' '34'}), 'unicode')  native2unicode(hex2dec({'00' 'a0'}), 'unicode') ', x ' native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'91' 'cf'}), 'unicode') ] [native2unicode(hex2dec({'90' '1f'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode')  native2unicode(hex2dec({'58' '34'}), 'unicode')  native2unicode(hex2dec({'00' 'a0'}), 'unicode') ', y ' native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'91' 'cf'}), 'unicode') ] [native2unicode(hex2dec({'90' '1f'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode')  native2unicode(hex2dec({'58' '34'}), 'unicode')  native2unicode(hex2dec({'00' 'a0'}), 'unicode') ', z ' native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'91' 'cf'}), 'unicode') ]});

model.result.export('data1').set('filename', s_model_path + "/" + s_mesh_nodes_filename);
