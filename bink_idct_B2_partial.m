function out=bink_idct_B2_partial(in,stages)
  % extract rows (with input permutation)
  c0 = in(1,:);
  d4 = in(2,:);
  c2 = in(3,:);
  d6 = in(4,:);
  c1 = in(5,:);
  d5 = in(6,:);
  c3 = in(7,:);
  d7 = in(8,:);

  % odd stage 4
  c4 = d4;
  c5 = d5 + d6;
  c7 = d5 - d6;
  c6 = d7;

  if stages > 1
    % odd stage 3
    b4 = c4 + c5;
    b5 = c4 - c5;
    b6 = c6 + c7;
    b7 = c6 - c7;

    % even stage 3
    b0 = c0 + c1;
    b1 = c0 - c1;
    b2 = c2 + c2/4 + c3/2;
    b3 = c2/2 - c3 - c3/4;

    if stages > 2
      % odd stage 2
      a4 = b7/4 + b4 + b4/4 - b4/16;
      a7 = b4/4 - b7 - b7/4 + b7/16;
      a5 = b5 - b6 + b6/4 + b6/16;
      a6 = b6 + b5 - b5/4 - b5/16;

      % even stage 2
      a0 = b0 + b2;
      a1 = b1 + b3;
      a2 = b1 - b3;
      a3 = b0 - b2;

      if stages > 3
        % stage 1
        o0 = a0 + a4;
        o1 = a1 + a5;
        o2 = a2 + a6;
        o3 = a3 + a7;
        o4 = a3 - a7;
        o5 = a2 - a6;
        o6 = a1 - a5;
        o7 = a0 - a4;
      else
        o0 = b0;
        o1 = b1;
        o2 = b2;
        o3 = b3;
        o4 = b4;
        o5 = b5;
        o6 = b6;
        o7 = b7;
      end
    else
      o0 = b0;
      o1 = b1;
      o2 = b2;
      o3 = b3;
      o4 = b4;
      o5 = b5;
      o6 = b6;
      o7 = b7;
    end
  else
    o0 = c0;
    o1 = c1;
    o2 = c2;
    o3 = c3;
    o4 = c4;
    o5 = c5;
    o6 = c6;
    o7 = c7;
  end

  % output
  out = [o0; o1; o2; o3; o4; o5; o6; o7];
end
