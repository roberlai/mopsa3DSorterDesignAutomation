function ...
  [r, y_init, y_init_shift, y_gap, x_init, x_gap, theta, phi, fixed_side] ...
    =  readRecipe(filename)

  fprintf(1,'%s\n',filename);

  fid = fopen(filename);
  title = "";
  title = fgetl(fid);

  vals = fscanf(fid, "%f %f %f %f %f %f %f %f %f\n");

  r = vals(1);

  theta = vals(2);
  phi = vals(3);

  y_init = vals(4);
  y_init_shift = vals(5);
  y_gap = vals(6);

  x_init = vals(7);
  x_gap = vals(8);

  fixed_side = vals(9);
end
