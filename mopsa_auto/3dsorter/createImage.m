
model.component('comp1').view('view1').set('transparency', true);
model.result.export.create('img1', 'Image');
model.result.export('img1').set('sourceobject', 'pg1');
model.result.export('img1').set('pngfilename', s_output_image);
model.result.export('img1').set('size', 'current');
model.result.export('img1').set('unit', 'px');
model.result.export('img1').set('height', '540');
model.result.export('img1').set('width', '885');
model.result.export('img1').set('lockratio', 'off');
model.result.export('img1').set('resolution', '96');
model.result.export('img1').set('antialias', 'on');
model.result.export('img1').set('zoomextents', 'on');
model.result.export('img1').set('fontsize', '9');
model.result.export('img1').set('customcolor', [1 1 1]);
model.result.export('img1').set('background', 'current');
model.result.export('img1').set('gltfincludelines', 'on');
model.result.export('img1').set('title1d', 'on');
model.result.export('img1').set('legend1d', 'on');
model.result.export('img1').set('logo1d', 'on');
model.result.export('img1').set('options1d', 'on');
model.result.export('img1').set('title2d', 'on');
model.result.export('img1').set('legend2d', 'on');
model.result.export('img1').set('logo2d', 'on');
model.result.export('img1').set('options2d', 'off');
model.result.export('img1').set('title3d', 'on');
model.result.export('img1').set('legend3d', 'on');
model.result.export('img1').set('logo3d', 'on');
model.result.export('img1').set('options3d', 'on');
model.result.export('img1').set('axisorientation', 'on');
model.result.export('img1').set('grid', 'on');
model.result.export('img1').set('axes1d', 'on');
model.result.export('img1').set('axes2d', 'on');
model.result.export('img1').set('showgrid', 'on');
model.result.export('img1').set('target', 'file');
model.result.export('img1').set('qualitylevel', '92');
model.result.export('img1').set('qualityactive', 'off');
model.result.export('img1').set('imagetype', 'png');
model.result.export('img1').set('lockview', 'off');
model.result.export('img1').run();

