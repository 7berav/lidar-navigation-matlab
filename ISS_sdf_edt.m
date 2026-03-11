%% 0) 준비
syms x y z
order_Data  = 4;
TermsC_Data = homogeneFischerTerms(order_Data);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
range_Data  = [-40 40 -35 35 -30 20];
%range_Data  = [-35 35 -10 18 -30 22];

faceColor_Data = [0.95, 0.72, 0.5];
faceAlpha_Data = 0.50; meshDen_Data = 100;
h = 0.1;    

xmin = range_Data(1); xmax = range_Data(2);
ymin = range_Data(3); ymax = range_Data(4);
zmin = range_Data(5); zmax = range_Data(6);

x = xmin:h:xmax;  y = ymin:h:ymax;  z = zmin:h:zmax;
[Nx, Ny, Nz] = deal(numel(x), numel(y), numel(z));

meta.origin  = [xmin, ymin, zmin];
meta.spacing = [h, h, h];
meta.size    = [Nx, Ny, Nz];

%% ===================== 1) 데이터 로드 =====================
if exist('DATA4.mat','file')
    S_Data = load('DATA4.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
else
    error('DATA4.mat 파일이 없습니다.');
end

hasFields = all(isfield(Data_Data, {'Beta','Disp'}));
if ~hasFields, error('Data_Data에 Beta/Disp 필드가 없습니다.'); end
valid = arrayfun(@(s) ~isempty(s.Beta) && ~isempty(s.Disp), Data_Data);
idxAll = find(valid); 
if isempty(idxAll), error('유효 항목 없음'); end

outDir = 'func_bank'; if ~exist(outDir,'dir'), mkdir(outDir); end
% m 파일들이 있는 폴더를 경로에 추가
if ~exist('outDir','var'), outDir = 'func_bank'; end
if exist(outDir,'dir'), addpath(outDir); end
if ~exist(outDir,'dir')
    warning('outDir가 없습니다: %s (이 인덱스는 스킵됩니다)', outDir);
    % 여기서 return 하거나, 루프에서 continue 하세요.
end
% 사용할 도형 인덱스(예: 1:25 전체 또는 부분집합)
shape_ids = idxAll(:)';       % 필요에 맞게 수정 (예: 1:25)

% 각 도형의 함수핸들 미리 로드(파라포 대비)
G = cell(1, numel(shape_ids));      % g 핸들
for t = 1:numel(shape_ids)
    j = shape_ids(t);
    mpath = sprintf('f_ISS_%d.mat', j);
    if ~exist(mpath,'file')
        warning('mat파일을 찾지 못했습니다: %s (스킵)', mpath);
        G{t} = [];
        continue;
    end
    S = load(mpath, 'f_shifted_handle');
    if ~isfield(S,'f_shifted_handle') || ~isa(S.f_shifted_handle,'function_handle')
        warning('i=%d: %s 에 f_shifted_handle 함수핸들이 없습니다. (스킵)', j, mpath);
        G{t} = [];
        continue;
    end
    G{t} = S.f_shifted_handle;   % g(X,Y,Z)
end

% 비어있는 것 제거
keep = ~cellfun(@isempty, G);
G = G(keep);
shape_ids = shape_ids(keep);
nShape = numel(shape_ids);
if nShape == 0, error('유효 g 핸들이 하나도 없습니다.'); end
fprintf('유효 도형 수: %d\n', nShape);

%% ===================== 2) 내부 마스크(합집합) 및 내부 라벨 구축 =====================
% 큰 배열: 논리 inside, 내부 라벨(겹침 시 |φ| 최소 도형)
inside = false(Nx,Ny,Nz);                 % 논리
label  = zeros(Nx,Ny,Nz, 'uint8');       % 내부 보셀의 "가장 가까운" 도형 id

ncores = feature('numcores');
p = gcp('nocreate');
if isempty(p)
    p = parpool('local', ncores);   % 또는 parpool('local', max(1,ncores-1))
end

Gc = parallel.pool.Constant(G); 
ids = shape_ids; 
% 병렬화: z-슬라이스 병렬
parfor k = 1:Nz
    Zk = z(k);
    [XX,YY] = ndgrid(x,y);        % (Nx, Ny)
    inside_k = false(Nx,Ny);      % 로컬
    label_k  = zeros(Nx,Ny,'uint16');
    minAbsPhi_k = inf(Nx,Ny,'single');
    Glocal = Gc.Value;
    % 모든 도형을 순회하며 φ 평가 → 내부 갱신
    for t = 1:nShape
        g = Glocal{t};
        if isempty(g), continue; end
        % φ = g - 1
        phi = single(g(XX,YY, Zk) - 1);   % g가 벡터화되어 있다고 가정
        mask = (phi <= 0);
        if any(mask(:))
            inside_k = inside_k | mask;
            % 내부 겹침에서 |φ|가 더 작은 도형으로 라벨 갱신
            absphi = abs(phi);
            better = mask & (absphi < minAbsPhi_k);
            if any(better(:))
                minAbsPhi_k(better) = absphi(better);
                label_k(better) = uint16(ids(t));
            end
        end
    end

    inside(:,:,k) = inside_k;
    label(:,:,k)  = label_k;   % 내부가 아닌 곳은 0
end

%% ===================== 3) 거리변환(EDT) → 부호있는 거리장 =====================
% 외부(+): 외부 보셀에서 내부까지 거리
% 내부(-): 내부 보셀에서 외부까지 거리
% bwdist는 격자 단위 거리를 반환하므로 마지막에 h를 곱함.
[Dout, idxOut] = bwdist( inside, 'euclidean' );
[Din,  idxIn ] = bwdist(~inside, 'euclidean' );

SDF = single(Dout);                % 외부는 +
SDF(inside) = -single(Din(inside));
SDF = SDF * h;                     % 물리 단위(해상도 h) 반영

nearestId = zeros(Nx,Ny,Nz, 'uint16');
% 외부: idxOut이 가리키는 위치의 내부 라벨을 복사
linearTrue = find(inside);                 % true 위치들의 선형 인덱스
Lab = zeros(Nx*Ny*Nz,1,'uint16');          % 라벨의 선형화 복제본
Lab(linearTrue) = reshape(label(inside), [], 1);   % true 위치만 보유
nearestId(~inside) = reshape( Lab(idxOut(~inside)), [], 1 );
% 내부: 그 보셀의 label (없으면 0)
nearestId(inside) = label(inside);
%%
% 관심 z 인덱스 선택 (중앙 슬라이스 예시)
k0 = round(numel(z)/2);
k0 = 340;
figure; 
[X,Y,Z] = ndgrid(x,y,z);

% (a) Dout 슬라이스 (외부에서 의미 있음)
imagesc(X,Y, squeeze(Dout(:,:,k0))' * h); axis image; colorbar;
title(sprintf('Dout @ z=%.3f (물리단위)', z(k0)));

fv = isosurface(X,Y,Z,SDF, 0);  % 0-레벨셋
figure; p = patch(fv); 
p.FaceColor = [0.9 0.7 0.5]; p.EdgeColor = 'none'; p.FaceAlpha = 0.3;
axis equal; camlight headlight; lighting gouraud; view(3);
title('SDF = 0 (표면)');
xlabel('X(cm)');ylabel('Y (cm)');zlabel('Z (cm)');
% 임의 등거리면 (예: +1, +5)
hold on;
isosurface(X,Y,Z,SDF, +0.5); %isosurface(SDF, -0.2);
hold off;
%% ===================== 4) 저장 =====================
outDir = 'sdf_bank';
if ~exist(outDir,'dir'), mkdir(outDir); end
outPath = fullfile(outDir, 'sdf_union_ISS.mat');
save(outPath, 'SDF', 'nearestId', 'x', 'y', 'z', 'meta', '-v7.3');
fprintf('저장 완료: %s\n', outPath);