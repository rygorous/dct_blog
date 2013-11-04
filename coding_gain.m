function g = coding_gain(fwd_xform, rho)
  % autocorrelation matrix of first-order gauss-markov process
  % with autocorrellation coefficient rho and unit variance
  len = size(fwd_xform, 1);
  R = toeplitz(rho .^ (0:len-1));

  % variance of transform coefficients
  sq_var_xform = diag(fwd_xform * R * fwd_xform');

  % squared norms of rows of inv_xform
  inv_xform = inv(fwd_xform);
  sq_norms = diag(inv_xform' * inv_xform);
  gain = 1.0 / mean(sq_var_xform .* sq_norms, 'g');
  g = 10*log10(gain);
end
