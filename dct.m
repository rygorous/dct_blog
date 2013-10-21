function r = dct(in)
    % computes the DCT-II of the input vector(s)
    % (this is the version with orthogonal normalization)
    N = size(in,1);
    r = diag([sqrt(1/N) sqrt(2/N)*ones(1,N-1)]) * cos(pi/N * (0:N-1)' * ((0:N-1) + 0.5)) * in;
end
