function out=dct_vc1(in)
  % source: multimedia wiki
  % http://wiki.multimedia.cx/index.php?title=VC-1_Tables#Transform_Matrix_8x8
  M = [
    12,  12,  12,  12,  12,  12,  12,  12;
    16,  15,   9,   4,  -4,  -9, -15, -16;
    16,   6,  -6, -16, -16,  -6,   6,  16;
    15,  -4, -16,  -9,   9,  16,   4, -15;
    12, -12, -12,  12,  12, -12, -12,  12;
     9, -16,   4,  15, -15,  -4,  16,  -9;
     6, -16,  16,  -6,  -6,  16, -16,   6;
     4,  -9,  15, -16,  16, -15,   9,  -4;
  ];

  out = M * in;
end
