function out=dct_h264(in)
  I4 = eye(4);
  J4 = fliplr(I4);

  % pass 1 butterflies
  pass1 = [I4 J4;J4 -I4];

  % even pass 2
  I2 = eye(2);
  J2 = fliplr(I2);
  even2 = [I2 J2;J2 -I2];
  
  % odd pass 2
  odd2 = 0.5*[0 2 2 3;-2 -3 0 2;2 0 -3 2;3 -2 2 0];

  pass2 = blkdiag(even2, odd2);

  % even-even pass 3
  eveneven3 = [1 1;1 -1];

  % even-odd pass 3
  evenodd3 = 0.5*[1 2;-2 1];

  % odd pass 3
  odd3 = 0.25*[4 0 0 1;0 4 1 0;0 -1 4 0;1 0 0 -4];

  pass3 = blkdiag(eveneven3, evenodd3, odd3);

  % output permute
  permute = eye(8) ([1 5 3 6 2 7 4 8],:);

  out = permute * pass3 * pass2 * pass1 * in;
end

% diag(M*M') = [8, 9+1/32, 5, 9+1/32, 8, 9+1/32, 5, 9+1/32]
