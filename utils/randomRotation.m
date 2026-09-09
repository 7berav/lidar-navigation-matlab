function R = randomRotation()
%RANDOMROTATION  Uniformly random 3-D rotation matrix (Rodrigues form).
%   R = randomRotation() draws a random unit axis and a random angle in
%   [0, 2*pi) and returns the corresponding 3-by-3 rotation matrix.
%
%   Draws from the ambient RNG, so seeding with rng(...) before the call
%   makes the view reproducible -- Sec 3.1 and Sec 3.2 rely on this to pair
%   the same shape realisations across sections.
%
%   Shared by Sec 3.2 scripts -- do not fork.

    a = randn(3,1); a = a / norm(a);
    ang = 2*pi*rand;
    K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    R = eye(3) + sin(ang)*K + (1-cos(ang))*(K*K);
end
