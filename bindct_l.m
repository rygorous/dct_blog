function out=bindct_l(in)
  I4 = eye(4);
  J4 = fliplr(I4);
  I2 = eye(2);
  J2 = fliplr(I2);

  % full-precision values
  p1 = 0.4142135623;
  u1 = 0.3535533905;
  p2 = 0.3033466836;
  u2 = 0.5555702330;
  p3 = 0.3033466836;
  p4 = 0.0984914033;
  u3 = 0.1950903220;
  p5 = 0.0984914033;

  % binDCT-L1
  p1 = 13/32;
  u1 = 11/32;
  p2 = 19/64;
  u2 = 9/16;
  p3 = 19/64;
  p4 = 3/32;
  u3 = 3/16;
  p5 = 3/32;

  % binDCT-L5
  %p1 = 1/2;
  %u1 = 1/2;
  %p2 = 1/4;
  %u2 = 1/2;
  %p3 = 1/4;
  %p4 = 1/8;
  %u3 = 1/4;
  %p5 = 1/8;

  % stage 1
  stage1 = [I4 J4;J4 -I4];

  % stage 2
  stage2even = [I2 J2;J2 -I2];

  oddrot1 = [1 0;-p3 1] * [1 u2;0 1] * [1 0;-p2 1];
  oddrot2 = [1 0;-p5 1] * [1 u3;0 1] * [1 0;-p4 1];

  perm_in = eye(4)([1 4 2 3],:);
  stage2odd = perm_in' * blkdiag(oddrot1, oddrot2) * perm_in;

  stage2 = blkdiag(stage2even, stage2odd);

  % stage 3
  stage3eveneven = [1 0;1/2 -1] * [1 1;0 1];
  stage3evenodd = [1 0;-u1 1] * [-1 p1;0 1];

  stage3odd = [1 0 1 0;0 -1 0 1;1 0 -1 0;0 1 0 1];
  stage3 = blkdiag(stage3eveneven, stage3evenodd, stage3odd);

  % stage 4
  stage4odd = perm_in' * blkdiag([-1 1/2;0 1] * [1 0;1 1], I2) * perm_in;
  stage4 = blkdiag(I4, stage4odd);

  permute = eye(8) (:,([0 4 6 2 7 3 5 1] + 1));
  % normalization (if desired)
  % permute = permute * diag([sin(pi/4)/2 sin(pi/4) sin(3*pi/8)/2 1/(2*sin(3*pi/8)) 1/sqrt(2) 1/2 1/2 1/sqrt(8)]);

  out = permute * stage4 * stage3 * stage2 * stage1 * in;
end
