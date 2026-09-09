%% report_tables.m — 3.2 결과 CSV 를 표로 출력 (재실행 없음, 읽기 전용)
%
% 스윕을 다시 돌리지 않고 최신(또는 지정) CSV 에서 필요한 슬라이스를 뽑는다.
%
% 사용:
%   run('experiments/sec3_2_convergence/report_tables.m')          % 최신 eval fold
%   useFold='tuning'; run(...)                                     % λ 선택용 fold
%   csvPath='test_ridge/sec32_single_..../ransac_outlier_results.csv'; run(...)
%   pickShape='cube'; pickOcc=0.5; run(...)                        % 슬라이스 변경
%
% 워크스페이스 override: csvPath, useFold, pickShape, pickEps, pickOcc, useQuick
%
% 표 읽기 규약:
%   theta_axis : 각도(도). 무작위 추정의 중앙값이 60도 = 정보 없음 기준선.
%                전체 시드로 집계한다(실패 제외 시 생존편향).
%   그 외 지표 : 실패 시드 제외 + failRate 별도 보고.
%   cd_occ     : occ_empty 와 같이 볼 것. occ_empty 가 높으면 실거리가 아니라
%                bbox 대각 페널티다.

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))

if ~exist('useQuick','var')||isempty(useQuick), useQuick=false; end
repoRoot_=fullfile(fileparts(mfilename('fullpath')),'..','..');
% run() 이 스크립트 폴더로 cd 하므로 상대경로는 저장소 루트 기준으로 다시 푼다.
if exist('csvPath','var') && ~isempty(csvPath) && ~isfile(csvPath)
    if isfile(fullfile(repoRoot_,csvPath)), csvPath=fullfile(repoRoot_,csvPath); end
end
% useFold: 'eval'(기본, 논문 수치 = seedOffset 100) / 'tuning'(λ 선택용) / 'latest'
% 주의: 'latest' 는 가장 최근 런을 집으므로 fold 가 섞일 수 있다. 보고용은 'eval'.
if ~exist('useFold','var')||isempty(useFold), useFold='eval'; end
if ~exist('csvPath','var')||isempty(csvPath)
    outRoot=fullfile(repoRoot_,'test_ridge');
    dd=dir(outRoot); nm={dd([dd.isdir]).name};
    hits=nm(startsWith(nm,'sec32_single_'));
    if ~useQuick, hits=hits(~contains(hits,'_quick')); end
    % 진행 중인 런은 폴더만 있고 CSV 가 없다. CSV 있는 것만 후보로 둔다.
    hits=hits(cellfun(@(h) isfile(fullfile(outRoot,h,'ransac_outlier_results.csv')),hits));
    switch lower(useFold)
        case 'eval',   pick=hits(contains(hits,'_fold'));
        case 'tuning', pick=hits(~contains(hits,'_fold'));
        otherwise,     pick=hits;
    end
    if isempty(pick)
        warning('report_tables:noFold','%s fold 결과가 없어 최신 런을 사용합니다.',useFold);
        pick=hits;
    end
    assert(~isempty(pick),'sec32_single_* 결과 폴더가 없습니다.');
    pick=sort(pick); csvPath=fullfile(outRoot,pick{end},'ransac_outlier_results.csv');
end
T=readtable(csvPath);
fprintf('\n소스: %s\n행 %d개, 시드 %s\n', csvPath, height(T), mat2str([min(T.seed) max(T.seed)]));

% 구버전 CSV 호환
has=@(c) ismember(c,T.Properties.VariableNames);
if ~has('occ'),        T.occ        = 0.30*ones(height(T),1); end
if ~has('occ_empty'),  T.occ_empty  = nan(height(T),1);       end
if ~has('otype'),      T.otype      = repmat({'uniform'},height(T),1); end

% 실패 기준: cd_all 문턱(스크립트에서 이미 반영) 또는 중심오차 15%*bbox
T.failX = double(T.fail==1 | hypot(T.e_par,T.e_perp) > 0.15);

if ~exist('pickShape','var')||isempty(pickShape), pickShape='hexagon'; end
if ~exist('pickEps','var')  ||isempty(pickEps),   pickEps=0.2;         end
if ~exist('pickOcc','var')  ||isempty(pickOcc),   pickOcc=0.30;        end
mIt=mode(T.maxIter);

sel=@(s,e,l,c,o) strcmp(T.shape,s)&abs(T.eps-e)<1e-9&abs(T.lambda-l)<1e-12 & ...
                 abs(T.occ-c)<1e-9&T.maxIter==mIt&strcmp(T.otype,o);
med  =@(m,v) median(T.(v)(m&T.failX==0),'omitnan');   % 실패 제외
medA =@(m,v) median(T.(v)(m),'omitnan');              % 전체 시드
rate =@(m,v) 100*mean(T.(v)(m),'omitnan');

lams=unique(T.lambda)'; occs=unique(T.occ)'; epss=unique(T.eps)'; shps=unique(T.shape)';

%% ---- [1] λ × 도형 (occ=pickOcc, eps=pickEps) ----
fprintf('\n===== [1] λ × 도형   (occ=%.2f, eps=%.2f, uniform) =====\n',pickOcc,pickEps);
for v={'F1','cd_all','cd_occ'}
    fprintf('\n[%s]  %-9s',v{1},'shape'); fprintf('%9s',arrayfun(@(l) string(lamName(l)),lams));
    fprintf('   %s\n',ternS(strcmp(v{1},'F1'),'argmax','argmin'));
    for si=1:numel(shps)
        y=arrayfun(@(l) med(sel(shps{si},pickEps,l,pickOcc,'uniform'),v{1}),lams);
        if all(isnan(y)), continue; end
        if strcmp(v{1},'F1'), [~,ix]=max(y); else, [~,ix]=min(y); end
        fprintf('      %-9s',shps{si}); fprintf('%9.3f',y); fprintf('   %s\n',lamName(lams(ix)));
    end
end

%% ---- [2] occ × λ  (eps=pickEps, pickShape) ----
fprintf('\n===== [2] occ × λ   (%s, eps=%.2f, uniform) =====\n',pickShape,pickEps);
for v={'cd_all','cd_occ','F1','failX','occ_empty','theta_axis'}
    lamHere=lams(arrayfun(@(l) sum(arrayfun(@(c) any(sel(pickShape,pickEps,l,c,'uniform')),occs))>=2,lams));
    if isempty(lamHere), continue; end
    fprintf('\n[%s]%s\n%-6s',v{1},ternS(strcmp(v{1},'theta_axis'),'  (도, 전체시드. 무작위=60도)', ...
        ternS(any(strcmp(v{1},{'failX','occ_empty'})),'  (%)','  (실패 제외 중앙값)')),'occ');
    fprintf('%10s',arrayfun(@(l) string(lamName(l)),lamHere)); fprintf('\n');
    for c=occs
        fprintf('%-6.2f',c);
        for l=lamHere
            m=sel(pickShape,pickEps,l,c,'uniform');
            if ~any(m), fprintf('%10s','-'); continue; end
            switch v{1}
                case {'failX','occ_empty'}, fprintf('%10.0f',rate(m,v{1}));
                case 'theta_axis',          fprintf('%10.1f',medA(m,v{1}));
                otherwise,                  fprintf('%10.3f',med(m,v{1}));
            end
        end
        fprintf('\n');
    end
end

%% ---- [3] occ × ε  (오염이 곡선을 평행이동시키는가 기울기를 바꾸는가) ----
epsHere=epss(arrayfun(@(e) sum(arrayfun(@(c) any(sel(pickShape,e,0,c,'uniform')),occs))>=2,epss));
if numel(epsHere)>1
    fprintf('\n===== [3] occ × ε   (%s, uniform) =====\n',pickShape);
    for v={'cd_all','failX'}
        fprintf('\n[%s]\n%-6s',v{1},'occ');
        for e=epsHere, for l=unique([0 1e-2]), fprintf('%14s',sprintf('e%.0f%% %s',100*e,lamName(l))); end, end
        fprintf('\n');
        for c=occs
            fprintf('%-6.2f',c);
            for e=epsHere
                for l=unique([0 1e-2])
                    m=sel(pickShape,e,l,c,'uniform');
                    if ~any(m), fprintf('%14s','-');
                    elseif strcmp(v{1},'failX'), fprintf('%14.0f',rate(m,v{1}));
                    else, fprintf('%14.3f',med(m,v{1})); end
                end
            end
            fprintf('\n');
        end
    end
end

%% ---- [4] 손실 분해 사다리 (eps 축) ----
if has('cd_occ_oracle')
    fprintf('\n===== [4] 손실 분해 (%s, occ=%.2f, λ=1e-2, uniform) =====\n',pickShape,pickOcc);
    fprintf('%-6s %11s %11s %11s | %10s %10s\n','eps','L2d direct','L2 oracle','L3 RANSAC','L3-L2','L2d-L2');
    for e=epss
        m=sel(pickShape,e,1e-2,pickOcc,'uniform'); if ~any(m), continue; end
        d=med(m,'cd_occ_direct'); o=med(m,'cd_occ_oracle'); r=med(m,'cd_occ');
        fprintf('%-6.2f %11.4f %11.4f %11.4f | %+10.4f %+10.4f\n',e,d,o,r,r-o,d-o);
    end
end

%% ---- [5] 자세: 적합 vs raw PCA ----
fprintf('\n===== [5] 장축 자세 (occ=%.2f, eps=%.2f, λ=1e-2) — 무작위 추정 중앙값 = 60도 =====\n',pickOcc,pickEps);
fprintf('%-9s %9s %9s %11s %11s %s\n','shape','gap GT','gap fit','theta fit','theta PCA','유효시드');
for si=1:numel(shps)
    m=sel(shps{si},pickEps,1e-2,pickOcc,'uniform'); if ~any(m), continue; end
    fprintf('%-9s %9.2f %9.2f %11s %11s %d/%d\n',shps{si}, ...
        medA(m,'axis_gap_gt'), medA(m,'axis_gap_fit'), ...
        num2str(medA(m,'theta_axis'),'%.1f'), num2str(medA(m,'theta_axis_pca'),'%.1f'), ...
        sum(~isnan(T.theta_axis(m))), sum(m));
end
fprintf('\n(gap<=1.15 이면 2차 모멘트가 등방이라 자세 평가 불가 → NaN. cube/sphere 가 여기 해당.)\n');

%% ---- [6] 오염 유형 ----
otys=unique(T.otype)';
if numel(otys)>1
    fprintf('\n===== [6] 오염 유형 × ε  (%s, occ=%.2f, λ=1e-2) =====\n',pickShape,pickOcc);
    fprintf('%-6s',' eps'); for o=otys, fprintf('%26s',sprintf('%s  F1/prec/CD/fail%%',o{1})); end; fprintf('\n');
    for e=epss
        fprintf('%-6.2f',e);
        for o=otys
            m=sel(pickShape,e,1e-2,pickOcc,o{1});
            if ~any(m), fprintf('%26s','-'); continue; end
            fprintf('%26s',sprintf('%.2f/%.2f/%.3f/%.0f',med(m,'F1'),med(m,'precision'),med(m,'cd_all'),rate(m,'failX')));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% ============================ local functions ============================
function s=lamName(l), if l==0, s='OLS'; else, s=sprintf('%.0e',l); end, end
function s=ternS(c,a,b), if c, s=a; else, s=b; end, end
