function [phiFun, gradStackFun, nT] = makePuncsGradStack(Terms)
    % Terms: sym(1 x nT) — 사용자가 고른 항
    syms x y z
    nT   = numel(Terms);
    phi  = Terms(:).';                    % 1 x nT
    dpx  = diff(phi, x);                  % 1 x nT
    dpy  = diff(phi, y);                  % 1 x nT
    dpz  = diff(phi, z);                  % 1 x nT

    % N×3 입력(XYZ)에 대해 N×nT 반환
    phiFun = matlabFunction(phi, 'Vars',{[x y z]}, 'Optimize',true);

    % gradStackFun(XYZ) → (3N)×nT, 블록순서 [dx; dy; dz]
    dpxFun = matlabFunction(dpx, 'Vars',{[x y z]}, 'Optimize',true);
    dpyFun = matlabFunction(dpy, 'Vars',{[x y z]}, 'Optimize',true);
    dpzFun = matlabFunction(dpz, 'Vars',{[x y z]}, 'Optimize',true);

    gradStackFun = @(XYZ) [ dpxFun(XYZ); dpyFun(XYZ); dpzFun(XYZ) ];
end