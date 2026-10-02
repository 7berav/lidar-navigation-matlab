function [groupId, groupNames, partId, partNames, dist] = issPartLabels(Pq, glbFile, cacheFile)
% ISS 점군의 부품 정답 라벨 (평가 전용 — 분할 방법 안에서는 쓰지 않는다)
%
% ISS_stationary.xyz 는 ISS_stationary.glb 메시를 좌표계 그대로 내보낸 점군이다.
% GLB 의 삼각형 면을 면적 비례로 촘촘히 샘플링하고, 각 질의점에 가장 가까운 샘플의
% 부품(메시 노드)을 붙인다. 부품은 노드 이름 앞의 조립 번호("08 P6 Truss_01" → 08)로 묶는다.
%
% 입력
%   Pq        : M x 3 질의점 (buildGraph 를 거친 점군)
%   glbFile   : GLB 경로                         (기본 'ISS_stationary.glb')
%   cacheFile : 표면 샘플 캐시(.mat)             (기본 out_cutting_ISS/gt/glb_samples.mat)
%
% 출력
%   groupId    : M x 1 부품 그룹 번호
%   groupNames : 그룹 이름 (string)
%   partId     : M x 1 부품(메시 노드) 번호
%   partNames  : 부품 이름 (string)
%   dist       : M x 1 가장 가까운 표면 샘플까지의 거리

if nargin < 2 || isempty(glbFile),   glbFile = 'ISS_stationary.glb'; end
if nargin < 3 || isempty(cacheFile), cacheFile = fullfile('out_cutting_ISS','gt','glb_samples.mat'); end

persistent Mdl Spart pNames gOfPart gNames cacheKey
key = [glbFile '|' cacheFile];
if isempty(Mdl) || ~strcmp(cacheKey, key)
    if exist(cacheFile, 'file')
        S = load(cacheFile);
    else
        S = sampleGLB(glbFile);
        d = fileparts(cacheFile);
        if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
        save(cacheFile, '-struct', 'S', '-v7.3');
    end
    Mdl = createns(double(S.Sx)); Spart = double(S.Spart); pNames = S.partNames;
    [gOfPart, gNames] = groupParts(pNames);
    cacheKey = key;
end
[i, dist] = knnsearch(Mdl, Pq);
partId     = Spart(i);
groupId    = gOfPart(partId);
groupNames = gNames;
partNames  = pNames;
end

function S = sampleGLB(glbFile)
% GLB(바이너리 glTF 2.0) → 월드좌표 삼각형 → 면적 비례 표면 샘플
fid = fopen(glbFile, 'r', 'ieee-le');
assert(fid > 0, 'GLB 파일을 열 수 없습니다: %s', glbFile);
hdr = fread(fid, 3, 'uint32');
assert(hdr(1) == 1179937895, 'GLB 형식이 아닙니다');          % 'glTF'
c0  = fread(fid, 2, 'uint32'); js = fread(fid, c0(1), 'uint8=>char')';
c1  = fread(fid, 2, 'uint32'); bin = fread(fid, c1(1), 'uint8=>uint8');
fclose(fid);

j = jsondecode(js);
nodes = aslist(j.nodes); meshes = aslist(j.meshes);
acc = aslist(j.accessors); bvs = aslist(j.bufferViews);
nN = numel(nodes);
parent = zeros(nN,1);
for i = 1:nN
    ch = getf(nodes{i}, 'children', []);
    if ~isempty(ch), parent(ch + 1) = i; end
end

Vc = {}; Fc = {}; Tc = {}; names = strings(0,1); pid = 0; base = 0;
for i = 1:nN
    nd = nodes{i};
    mi = getf(nd, 'mesh', []);
    if isempty(mi), continue; end
    pid = pid + 1; names(pid,1) = string(getf(nd, 'name', ''));
    M = localMat(nd); q = i;
    while parent(q) > 0, q = parent(q); M = localMat(nodes{q}) * M; end
    prims = aslist(meshes{mi + 1}.primitives);
    for p = 1:numel(prims)
        pr = prims{p};
        if getf(pr, 'mode', 4) ~= 4, continue; end            % 삼각형만
        pos = readAcc(bin, acc, bvs, pr.attributes.POSITION);
        pw  = pos * M(1:3,1:3).' + M(1:3,4).';
        ii  = getf(pr, 'indices', []);
        if isempty(ii), idx = (0:size(pos,1)-1)'; else, idx = readAcc(bin, acc, bvs, ii); end
        f = reshape(idx + 1 + base, 3, []).';
        Vc{end+1} = pw; Fc{end+1} = f; Tc{end+1} = pid * ones(size(f,1),1); %#ok<AGROW>
        base = base + size(pw,1);
    end
end
V = vertcat(Vc{:}); F = vertcat(Fc{:}); TM = vertcat(Tc{:});

st = rng; rng(0);                                           % 호출자의 난수 상태는 보존
A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
area = 0.5 * vecnorm(cross(B - A, C - A, 2), 2, 2);
a0 = max(sum(area) / 3e6, 1e-3);                            % 샘플 1개당 면적
ns = max(1, ceil(area / a0));
t  = repelem((1:size(F,1))', ns);
r1 = sqrt(rand(numel(t),1)); r2 = rand(numel(t),1);
Sx = (1-r1).*A(t,:) + r1.*(1-r2).*B(t,:) + r1.*r2.*C(t,:);
rng(st);
S = struct('Sx', single(Sx), 'Spart', uint16(TM(t)), 'partNames', names);
end

function A = readAcc(bin, acc, bvs, ai)
a = acc{ai + 1}; bv = bvs{a.bufferView + 1};
switch a.type
    case 'SCALAR', nc = 1; case 'VEC2', nc = 2; case 'VEC3', nc = 3; case 'VEC4', nc = 4;
    otherwise, error('지원하지 않는 accessor type: %s', a.type);
end
switch a.componentType
    case 5126, cls = 'single'; sz = 4;
    case 5125, cls = 'uint32'; sz = 4;
    case 5123, cls = 'uint16'; sz = 2;
    case 5121, cls = 'uint8';  sz = 1;
    otherwise, error('지원하지 않는 componentType: %d', a.componentType);
end
start  = getf(bv, 'byteOffset', 0) + getf(a, 'byteOffset', 0);
stride = getf(bv, 'byteStride', 0);
n = a.count;
if stride == 0 || stride == nc*sz
    raw = bin(start+1 : start + n*nc*sz);
else
    ix  = start + (0:n-1)'*stride + (1:nc*sz);
    raw = reshape(bin(ix.'), [], 1);
end
A = reshape(double(typecast(raw, cls)), nc, n).';
end

function M = localMat(nd)
m = getf(nd, 'matrix', []);
if ~isempty(m)
    M = reshape(m, 4, 4);                                   % glTF 는 열 우선
    return
end
t = getf(nd, 'translation', [0;0;0]); q = getf(nd, 'rotation', [0;0;0;1]); s = getf(nd, 'scale', [1;1;1]);
x = q(1); y = q(2); z = q(3); w = q(4);
R = [1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w);
     2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w);
     2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)];
M = eye(4); M(1:3,1:3) = R .* s(:).'; M(1:3,4) = t(:);
end

function [gOfPart, gName] = groupParts(nm)
% 부품 → 그룹: 이름 앞 두 자리 조립 번호로 묶고, 번호가 없는 것은 규칙으로 배정
key = strings(numel(nm),1);
for i = 1:numel(nm)
    tok = regexp(nm(i), '^(\d\d) ', 'tokens', 'once');
    if ~isempty(tok),                                              key(i) = tok{1};
    elseif startsWith(nm(i),'Destiny'),                            key(i) = "09";
    elseif startsWith(nm(i),'RapidScat'),                          key(i) = "45";
    elseif endsWith(nm(i),'.001') && startsWith(nm(i),'panel'),    key(i) = "05";   % Zvezda 태양전지판
    elseif startsWith(nm(i),'panel') || startsWith(nm(i),'hinge'), key(i) = "hinged solar panels (소속 미상)";
    else,                                                          key(i) = nm(i);
    end
end
[gKeys, ~, gOfPart] = unique(key, 'stable');
gName = strings(numel(gKeys),1);
for g = 1:numel(gKeys)
    if isempty(regexp(gKeys(g), '^\d\d$', 'once')), gName(g) = gKeys(g); continue; end
    cand = nm(gOfPart == g);
    num  = cand(~cellfun(@isempty, regexp(cellstr(cand), '^\d\d ', 'once')));
    if ~isempty(num), cand = num; end
    [~, k] = min(strlength(cand)); gName(g) = cand(k);
end
end

function c = aslist(x)
if iscell(x), c = x(:); elseif isstruct(x), c = num2cell(x(:)); else, c = {}; end
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
