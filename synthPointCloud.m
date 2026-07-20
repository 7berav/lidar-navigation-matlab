function [P, labels_gt, shapes] = synthPointCloud(shapes, noiseSigma)
% 합성 점군 생성: 여러 도형을 배치하고 정답 라벨과 함께 반환
%
% 입력
%   shapes     : struct 배열. 필드:
%                .type  'box' | 'ellipsoid' | 'hexprism' | 'cylinder' | 'sphere'
%                .N     점 개수
%                .size  도형 반치수 [a b c] (box/ellipsoid), 스칼라면 등방 스케일
%                .R     3x3 회전행렬 (생략 시 eye(3))
%                .t     1x3 평행이동 (생략 시 [0 0 0])
%   noiseSigma : 가우시안 노이즈 표준편차 (등방, 생략 시 0)
%
% 출력
%   P          : M x 3 점군 (전체 도형 합침)
%   labels_gt  : M x 1 uint32 정답 라벨 (도형 인덱스)
%   shapes     : 입력 shapes에 기본값 채워서 반환
%
% 예) 판 두 개가 45도로 교차:
%   s(1) = struct('type','box','N',3000,'size',[10 2 0.5]);
%   s(2) = struct('type','box','N',3000,'size',[10 2 0.5], ...
%                 'R',axang2rotm_local([0 0 1 pi/4]),'t',[0 0 1]);
%   [P, gt] = synthPointCloud(s, 0.02);

if nargin < 2, noiseSigma = 0; end

P = zeros(0,3);
labels_gt = zeros(0,1,'uint32');

for i = 1:numel(shapes)
    s = shapes(i);
    if ~isfield(s,'R') || isempty(s.R), s.R = eye(3); end
    if ~isfield(s,'t') || isempty(s.t), s.t = [0 0 0]; end
    if ~isfield(s,'size') || isempty(s.size), s.size = [1 1 1]; end
    if isscalar(s.size), s.size = s.size * [1 1 1]; end
    shapes(i) = s;

    switch lower(s.type)
        case 'box'
            Pi = sampleBoxSurface(s.N, s.size);
        case 'ellipsoid'
            Pi = generateRandomPointsOnEllipsoid(s.N, s.size(1), s.size(2), s.size(3));
        case 'hexprism'
            Pi = generateRandomPointsOnHexagonPrism(s.N) .* s.size;
        case 'cylinder'
            Pi = generateRandomPointsOnCylinder(s.N) .* s.size;
        case 'sphere'
            Pi = generateRandomPointsOnSphere(s.N) .* s.size;
        otherwise
            error('synthPointCloud: 알 수 없는 type "%s"', s.type);
    end

    Pi = Pi * s.R.' + s.t;
    P = [P; Pi]; %#ok<AGROW>
    labels_gt = [labels_gt; repmat(uint32(i), size(Pi,1), 1)]; %#ok<AGROW>
end

if noiseSigma > 0
    P = P + noiseSigma * randn(size(P));
end
end

function P = sampleBoxSurface(N, halfExt)
% 직육면체 표면 균일 샘플링 (면적 가중 → aspect ratio가 커도 밀도 균일)
a = halfExt(1); b = halfExt(2); c = halfExt(3);
% 면쌍별 면적: x면 2*(2b*2c), y면 2*(2a*2c), z면 2*(2a*2b)
areas = [b*c, a*c, a*b];
prob  = areas / sum(areas);
edges = cumsum([0 prob]);

u = rand(N,1);
faceAxis = discretize(u, edges);          % 1:x면, 2:y면, 3:z면
sgn = sign(rand(N,1) - 0.5); sgn(sgn==0) = 1;

P = zeros(N,3);
for ax = 1:3
    m = (faceAxis == ax);
    n = nnz(m);
    if n == 0, continue; end
    others = setdiff(1:3, ax);
    Q = zeros(n,3);
    Q(:,ax) = sgn(m) * halfExt(ax);
    Q(:,others(1)) = (2*rand(n,1)-1) * halfExt(others(1));
    Q(:,others(2)) = (2*rand(n,1)-1) * halfExt(others(2));
    P(m,:) = Q;
end
end
