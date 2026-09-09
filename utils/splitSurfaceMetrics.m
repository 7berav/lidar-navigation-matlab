function S = splitSurfaceMetrics(beta, PhiGrid, GX, GY, GZ, center, nSurf, PcleanFull, visMask, zThr, hdPctl)
%SPLITSURFACEMETRICS  CD/HD of an implicit fit against clean GT, split by visibility.
%   S = splitSurfaceMetrics(beta, PhiGrid, GX,GY,GZ, center, nSurf, ...
%                           PcleanFull, visMask, zThr, hdPctl)
%
%   Samples the f=1 level set (meshSurfacePoints) and scores it against the
%   occlusion-free GT point set PcleanFull in three scopes. Follows the Sec 3.1
%   convention exactly (lcurve_demo.m): the hidden-scope metric compares hidden
%   GT points against ONLY the fitted-surface points that fall in the hidden
%   region, not against the whole surface.
%
%     cd_occ/hd_occ : PcleanFull(~visMask,:)  vs  V(V(:,3) <  zThr,:)   <- occlusion headline
%     cd_vis/hd_vis : PcleanFull( visMask,:)  vs  V(V(:,3) >= zThr,:)
%     cd_all/hd_all : PcleanFull              vs  V
%     osf           : |V in hidden| / |V|  -- fitted-surface mass pushed into the
%                     unobserved region. Detects spurious lobes without a
%                     separate volumetric IoU.
%
%   A degenerate fit (empty level set) yields emptyPenalty = GT bbox diagonal
%   for every distance, and osf = NaN. S.V returns the sampled surface points
%   so callers can reuse them (e.g. principalAxis) without re-meshing.
%
%   Shared by Sec 3.1 and Sec 3.2 -- do not fork.

    emptyPenalty = norm(max(PcleanFull,[],1) - min(PcleanFull,[],1));

    V = meshSurfacePoints(beta, PhiGrid, GX, GY, GZ, center, nSurf);
    S.V = V;

    if isempty(V)
        [S.cd_occ, S.hd_occ] = deal(emptyPenalty);
        [S.cd_vis, S.hd_vis] = deal(emptyPenalty);
        [S.cd_all, S.hd_all] = deal(emptyPenalty);
        S.osf = NaN;
        S.empty = true; S.occ_empty = true;
        return;
    end
    S.empty = false;

    Vocc = V(V(:,3) <  zThr, :);
    Vvis = V(V(:,3) >= zThr, :);
    S.osf = size(Vocc,1) / max(size(V,1), 1);

    % cd_occ 가 실거리가 아니라 emptyPenalty 인가. 복원면이 hidden 영역에 아예
    % 들어가지 않으면(occ 낮을수록 흔함) cd_occ = bbox 대각이 찍힌다.
    % 이 플래그 없이 cd_occ 로 λ를 고르면 **강한 정규화가 유리해 보이는 착시**가 생긴다:
    % λ↑ 이면 등위면이 구(Sobolev 널스페이스 = r^order)가 되어 hidden 영역을 반드시
    % 채우므로 페널티를 회피한다. 3.1 실측에서 hexagon occ=0.2·λ=0 은 70% 가 페널티였고
    % 그래서 "cd_occ 최소 @ λ=3e0" 이 나왔다 — 적합 품질이 아니라 회피의 결과.
    % → cd_occ 를 집계할 때는 반드시 occ_empty 비율을 함께 보고할 것.
    S.occ_empty = isempty(Vocc);

    [S.cd_all, S.hd_all] = symmetricMetrics(PcleanFull,            V,    hdPctl, emptyPenalty);
    [S.cd_occ, S.hd_occ] = symmetricMetrics(PcleanFull(~visMask,:), Vocc, hdPctl, emptyPenalty);
    [S.cd_vis, S.hd_vis] = symmetricMetrics(PcleanFull( visMask,:), Vvis, hdPctl, emptyPenalty);
end
