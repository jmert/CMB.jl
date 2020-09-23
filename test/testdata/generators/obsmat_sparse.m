function obsmat_sparse()
  R = sparse([ ...
        0.1, 0.2, 0.0, 0.0, 0.0, 0.0; ...
        0.0, 0.0, 0.5, 0.5, 0.0, 0.0; ...
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0; ...
      ]);
  pixels_right.type = 'dummy_pixels';
  pixels_right.index = int64(reshape(0:5, [], 1));
  pixels_right.sub.extra = int8(1);
  pixels_left = int64(reshape(1:4, [], 1));

  [fdir,fname,fext] = fileparts(mfilename('fullpath'));
  [pdir,~,~] = fileparts(fdir);

  savename = fullfile(pdir, 'obsmat_sparse.mat');
  save(savename, '-v7.3', 'R', 'pixels_right', 'pixels_left');
end
