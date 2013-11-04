function r = dct_llm_basic(in)
    % 8-point DCT-II using the Loeffler-Ligtenberg-Moschytz factorization
    % as given in the paper
    % (i.e. orthogonal DCT-II scaled by sqrt(8))

    c1 = cos(pi/16);
    s1 = sin(pi/16);
    c3 = cos(3*pi/16);
    s3 = sin(3*pi/16);
    c6 = cos(6*pi/16);
    s6 = sin(6*pi/16);

    % stage 1 (odd/even separation butterflies)
    I4 = eye(4);
    J4 = fliplr(I4);
    stage1 = [I4 J4;J4 -I4];

    % stage 2
    I2 = eye(2);
    J2 = fliplr(I2);
    stage2_even = [I2 J2;J2 -I2];
    stage2_odd = [c3 0 0 s3; 0 c1 s1 0; 0 -s1 c1 0; -s3 0 0 c3];
    stage2 = blkdiag(stage2_even, stage2_odd);

    % stage 3
    % NOTE: in the paper, even stage 3 is written as having a "c1" rotation;
    % in fact it should be c6.
    stage3_even = blkdiag([1 1;1 -1], sqrt(2)*[c6 s6;-s6 c6]);
    stage3_odd = [1 0 1 0; 0 -1 0 1; 1 0 -1 0; 0 1 0 1];
    stage3 = blkdiag(stage3_even, stage3_odd);

    % stage 4
    stage4_odd = [-1 0 0 1; 0 sqrt(2) 0 0; 0 0 sqrt(2) 0; 1 0 0 1];
    stage4 = blkdiag(I4, stage4_odd);

    % NOTE: as given in the paper, several rows actually correspond to
    % *negated* DCT basis functions. fix this here.
    %stage4([4],:) = -stage4([4],:);

    % output permutation
    stage4 = stage4([1 8 3 6 2 7 4 5],:);

    r = stage4 * stage3 * stage2 * stage1 * in;
end
