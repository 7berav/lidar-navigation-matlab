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
    % makeVecFun을 로컬 함수로 인라인 (sym dispatch 우회)
    dpxFun = buildGradFun(dpx, x, y, z);
    dpyFun = buildGradFun(dpy, x, y, z);
    dpzFun = buildGradFun(dpz, x, y, z);

    Grads = @(XYZ) [ dpxFun(XYZ); dpyFun(XYZ); dpzFun(XYZ)];
end

function f = buildGradFun(dp, x, y, z)
    % dp: sym 1×nT  — 로컬 함수이므로 sym dispatch 문제 없음
    n  = numel(dp);
    fs = cell(1, n);
    for k = 1:n
        dpk = simplify(dp(k));
        if isequal(dpk, sym(0))
            fs{k} = @(XYZ) zeros(size(XYZ,1), 1);
        else
            fs{k} = matlabFunction(dpk, 'Vars', {[x y z]});
        end
    end
    f = @(XYZ) cell2mat(cellfun(@(g) g(XYZ), fs, 'UniformOutput', false));
end
