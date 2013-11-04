function out=bink_dct_approx(in)
  I4 = eye(4);
  J4 = fliplr(I4);
  I2 = eye(2);
  J2 = fliplr(I2);

  % all coding gains here computed for rho=0.95

  % integer, constrained-norm approximations to the rotations found by findorth.cpp

  %rot2 = 1/16 * [17 -7;7 17]; % variant a
  rot2 = 1/4 * [5 -2;2 5]; % variant b

  %rot1 = 1/8 * [8 -1;1 8]; % with rot2 variant a: 8.7971dB, 32A 10S transform.
  %rot3 = 1/8 * [7 4;-4 7]; % with rot2 variant b: 8.7968dB, 30A 10S transform.

  rot1 = 1/16 * [19 -4;4 19]; % with rot2 a: 8.8253dB, 38A 12S transform.
  rot3 = 1/16 * [16 11;-11 16]; % with rot2 b: 8.8250dB, 36A 12S transform.

  %rot1 = 1/64 * [65 -13;13 65]; % with rot2 a: 8.8259dB, 42A 15S transform.
  %rot3 = 1/64 * [55 37;-37 55]; % with rot2 b: 8.8255dB, 40A 15S transform.

  % build the matrices in order
  T8_0 = [I4 J4;I4 -J4];
  T8_0(8,:) = -T8_0(8,:); % sign flip so we get rotations not rotation-reflections
  
  T4_0 = [I2 J2;I2 -J2];
  T4_0(4,:) = -T4_0(4,:); % sign flip so we get rotations not rotation-reflections
  perm_in = eye(4)([1 4 2 3],:);
  T4_1 = perm_in' * blkdiag(rot1, rot3) * perm_in;

  T8_0_1 = blkdiag(T4_0, T4_1);

  C_II_2 = [1 1;1 -1];
  C_IV_2 = rot2;
  T8_0_1_0_0 = blkdiag(C_II_2, C_IV_2, C_II_2, C_II_2);

  A4 = [1 0 0 0; 0 1 0 1; 0 1 0 -1; 0 0 1 0];
  A8_0_1 = blkdiag(I4, A4);

  B8 = eye(8) (([0 4 2 6 1 5 3 7] + 1),:);

  out = B8 * A8_0_1 * T8_0_1_0_0 * T8_0_1 * T8_0 * in;
end
