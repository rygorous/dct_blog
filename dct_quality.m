function out=dct_quality(in_xform)
  ref_dct = dct(eye(size(in_xform)));

  % in_xform can be a scaled transform; determine the corresponding normalized transform.

  % first figure out gain factors
  gain_factors = diag(in_xform * in_xform');

  % we want to distribute gain equally between forward and inverse transform
  normalizer = diag(1 ./ sqrt(gain_factors));
  normalized_xform = normalizer * in_xform;

  % compute approximation error in matrix 2-norm
  % and coding gains for rho=0.90, rho=0.95
  l2 = norm(normalized_xform - ref_dct, 2);
  cg095 = coding_gain(normalized_xform, 0.95);
  cg090 = coding_gain(normalized_xform, 0.90);

  out = [l2, cg095, cg090];
end
