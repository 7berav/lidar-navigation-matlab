function [Funcs, Grads, nT] = makeFuncsGradsStack(Terms)
    % Terms: sym(1 x nT) — 사용자가 고른 항
    syms x y z
    nT   = numel(Terms);
    phi  = Terms(:).';                    % 1 x nT
    dpx  = diff(phi, x);                  % 1 x nT
    dpy  = diff(phi, y);                  % 1 x nT
    dpz  = diff(phi, z);                  % 1 x nT

    % ★ 0 → N×1로 강제
    dpx = dpx + 0*x;
    dpy = dpy + 0*y;
    dpz = dpz + 0*z;

    % N×3 입력(XYZ)에 대해 N×nT 반환
    Funcs = matlabFunction(phi, 'Vars',{[x y z]});

    % gradStackFun(XYZ) → (3N)×nT, 블록순서 [dx; dy; dz]
    %dpxFun = matlabFunction(dpx, 'Vars',{[x y z]});
    %dpyFun = matlabFunction(dpy, 'Vars',{[x y z]});
    dpxFun = makeVecFun(dpx, {[x y z]});
    dpyFun = makeVecFun(dpy, {[x y z]});
    dpzFun = makeVecFun(dpz, {[x y z]});

    %Grads = @(XYZ) [ dpxFun(XYZ); dpyFun(XYZ); dpzFun(XYZ) ]; % 3*nT
    Grads = @(XYZ) [ dpxFun(XYZ); dpyFun(XYZ); dpzFun(XYZ)];
end
