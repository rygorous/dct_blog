function out=bink_dct_B2(in)
  % extract rows
  i0 = in(1,:);
  i1 = in(2,:);
  i2 = in(3,:);
  i3 = in(4,:);
  i4 = in(5,:);
  i5 = in(6,:);
  i6 = in(7,:);
  i7 = in(8,:);

  % stage 1 - 8A
  a0 = i0 + i7;
  a1 = i1 + i6;
  a2 = i2 + i5;
  a3 = i3 + i4;
  a4 = i0 - i7;
  a5 = i1 - i6;
  a6 = i2 - i5;
  a7 = i3 - i4;

  % even stage 2 - 4A
  b0 = a0 + a3;
  b1 = a1 + a2;
  b2 = a0 - a3;
  b3 = a1 - a2;

  % even stage 3 - 6A 4S
  c0 = b0 + b1;
  c1 = b0 - b1;
  c2 = b2 + b2/4 + b3/2;
  c3 = b2/2 - b3 - b3/4;

  % odd stage 2 - 12A 8S
  % NB a4/4 and a7/4 are each used twice, so this really is 8 shifts, not 10.
  b4 = a7/4 + a4 + a4/4 - a4/16;
  b7 = a4/4 - a7 - a7/4 + a7/16;
  b5 = a5 + a6 - a6/4 - a6/16;
  b6 = a6 - a5 + a5/4 + a5/16;

  % odd stage 3 - 4A
  c4 = b4 + b5;
  c5 = b4 - b5;
  c6 = b6 + b7;
  c7 = b6 - b7;

  % odd stage 4 - 2A
  d4 = c4;
  d5 = c5 + c7;
  d6 = c5 - c7;
  d7 = c6;

  % permute/output
  out = [c0; d4; c2; d6; c1; d5; c3; d7];

  % total: 36A 12S
end
