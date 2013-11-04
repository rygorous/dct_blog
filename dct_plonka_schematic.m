function out=dct_plonka_schematic(in)
    I4 = eye(4);
    J4 = fliplr(I4);
    I2 = eye(2);
    J2 = fliplr(I2);
  
    c5 = cos(5*pi/16); s5 = sin(5*pi/16);
    c6 = cos(6*pi/16); s6 = sin(6*pi/16);
    c7 = cos(7*pi/16); s7 = sin(7*pi/16);

    % build the matrices in order
    T8_0 = [I4 J4;J4 -I4];
    
    T4_0 = [I2 J2;J2 -I2]; % scaled by sqrt(2)
    T4_1 = [c7 0 0 s7;0 c5 s5 0;0 -s5 c5 0;-s7 0 0 c7];
    T8_0_1 = blkdiag(T4_0, T4_1);
  
    C_II_2 = [1 1;1 -1]; % scaled by sqrt(2)
    C_IV_2 = [c6 s6;-s6 c6] * sqrt(2);
    T8_0_1_0_0 = blkdiag(C_II_2, C_IV_2, C_II_2, C_II_2);
  
    A4 = [sqrt(2) 0 0 0; 0 1 1 0; 0 1 -1 0; 0 0 0 -sqrt(2)];
    A8_0_1 = blkdiag(I4, A4);
  
    out = A8_0_1 * T8_0_1_0_0 * T8_0_1 * T8_0 * in;
    out = out([0 4 2 5 1 6 3 7]+1, :);
end
