function out=dct_plonka_paper(in)
    I4 = eye(4);
    J4 = fliplr(I4);
    I2 = eye(2);
    J2 = fliplr(I2);
  
    c1 = cos(1*pi/16); s1 = sin(1*pi/16);
    c2 = cos(2*pi/16); s2 = sin(2*pi/16);
    c3 = cos(3*pi/16); s3 = sin(3*pi/16);

    % build the matrices in order
    T8_0 = [I4 J4;I4 -J4];
    
    T4_0 = [I2 J2;I2 -J2]; % scaled by sqrt(2)
    T4_1 = [c1 0 0 s1;0 c3 s3 0;0 -s3 c3 0;s1 0 0 -c1];
    T8_0_1 = blkdiag(T4_0, T4_1);
  
    C_II_2 = [1 1;1 -1]; % scaled by sqrt(2)
    C_IV_2 = [c2 s2;s2 -c2] * sqrt(2);
    T8_0_1_0_0 = blkdiag(C_II_2, C_IV_2, C_II_2, C_II_2);
  
    A4 = [sqrt(2) 0 0 0; 0 1 0 1; 0 1 0 -1; 0 0 sqrt(2) 0];
    A8_0_1 = blkdiag(I4, A4);
  
    B8 = eye(8) (([0 4 2 6 1 5 3 7] + 1),:);
  
    out = B8 * A8_0_1 * T8_0_1_0_0 * T8_0_1 * T8_0 * in;
end
