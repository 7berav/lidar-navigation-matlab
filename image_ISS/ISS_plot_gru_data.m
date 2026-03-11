
%% 1) ISS 형상 표시
% 고정 파라미터
dt       = 0.1;                 % [s]
iterPlot = 4;                   % 3D 오버레이에 사용할 반복 index
sampleStep = 60;                % 6초 간격(0.1s * 60)
issDir   = 'image_ISS';         % ISS 관련 코드/데이터 폴더
figScriptName = 'figure_ISS_2'; % 형상 표시 스크립트/함수
sdfBaseName   = 'SDF_01';       % SDF 파일(또는 함수) 베이스 이름
outDir  = 'func_bank';          % 결과 저장 폴더
if ~exist(outDir,'dir'), mkdir(outDir); end

rangeSDF = [-40 40 -35 35 -30 20]; % [xmin xmax ymin ymax zmin zmax]
h = 0.1;                           % grid spacing

[Niter, Nt, ~] = size(deputy_rtn);


% 1) ISS 형상 fig 불러오기
f = openfig('image_ISS/figure_ISS_2.fig','reuse');
%%
hold on
plot3( ...
    squeeze(deputy_rtn(iterPlot,idx,1)), ...
    squeeze(deputy_rtn(iterPlot,idx,2)), ...
    squeeze(deputy_rtn(iterPlot,idx,3)), ...
    'o-', 'LineWidth', 1);
title('ISS surface + deputy trajectory (iter=1)')
axis equal
xlim([-20 20]);ylim([-25 25]);zlim([-10 10]);
view([0 0 1])
hold off
% 경로 추가
if exist(issDir,'dir'), addpath(issDir); end
SDF = [];
if exist('sdf_bank/SDF_02.mat','file')
    tmp = load('sdf_bank/SDF_02.mat');
    fn = fieldnames(tmp);
    SDF = tmp.(fn{1});
else
    error('SDF_02.mat 또는 SDF_02.m 없음')
end
%% statistic
idx = 1:sampleStep:Nt;
t = (idx-1)*dt;

mins = zeros(numel(idx),1);
meds = zeros(numel(idx),1);
means = zeros(numel(idx),1);
maxs = zeros(numel(idx),1);

xmin=rangeSDF(1); xmax=rangeSDF(2);
ymin=rangeSDF(3); ymax=rangeSDF(4);
zmin=rangeSDF(5); zmax=rangeSDF(6);

xvec = xmin:h:xmax;   % 길이 Nx
yvec = ymin:h:ymax;   % 길이 Ny
zvec = zmin:h:zmax;   % 길이 Nz

F = griddedInterpolant({xvec, yvec, zvec}, SDF, 'linear', 'none'); 

useGrid = ~isempty(SDF);
if useGrid, [Nx,Ny,Nz] = size(SDF); end
for k = 1:numel(idx)
    X = squeeze(deputy_rtn(:,idx(k),1)); X = X(:);
    Y = squeeze(deputy_rtn(:,idx(k),2)); Y = Y(:);
    Z = squeeze(deputy_rtn(:,idx(k),3)); Z = Z(:);

    % 범위 마스크(그리드/함수 공통)
    ix = 1 + round((X - xmin)/h);
    iy = 1 + round((Y - ymin)/h);
    iz = 1 + round((Z - zmin)/h);
    in = X>=xmin & X<=xmax & Y>=ymin & Y<=ymax & Z>=zmin & Z<=zmax;

    phi = nan(Niter,1);

    if useGrid
        % 그리드 인덱스도 유효해야 실제 조회 가능
        in = in & ix>=1 & ix<=Nx & iy>=1 & iy<=Ny & iz>=1 & iz<=Nz;
        ii = find(in);
        for n = ii.'
            phi(n,1) = F(X(n), Y(n), Z(n));
        end

    end

    % 집계
    v = phi(~isnan(phi));
    if isempty(v)
        mins(k)=NaN; meds(k)=NaN; means(k)=NaN; maxs(k)=NaN; stds(k)=NaN;
    else
        mins(k)=min(v);  meds(k)=median(v);  means(k)=mean(v);
        maxs(k)=max(v);  stds(k)=std(v,0);
    end
    N_oob(k) = sum(~in);          % 범위 밖(조회 불가로 간주)
    N_nan(k) = sum(isnan(phi));   % 최종 NaN 총계(= 범위 밖 + 평가실패)
end

%%
hold on
Xc = squeeze(deputy_rtn(iterPlot, idx, 1));
Yc = squeeze(deputy_rtn(iterPlot, idx, 2));
Zc = squeeze(deputy_rtn(iterPlot, idx, 3));

phi_c = nan(size(Xc));
if useGrid
    ix = 1 + round((Xc - xmin)/h);
    iy = 1 + round((Yc - ymin)/h);
    iz = 1 + round((Zc - zmin)/h);
    in = ix>=1 & ix<=Nx & iy>=1 & iy<=Ny & iz>=1 & iz<=Nz;
    for n = find(in)'
        phi_c(n) = SDF(ix(n), iy(n), iz(n));
    end
end
plot3(Xc, Yc, Zc, 'r-', 'LineWidth', 1);   % 검은색 선 (궤적)
scatter3(Xc, Yc, Zc, 10, phi_c, 'filled'); % 각 점 색상 = SDF 값
colormap(jet); colorbar;
%% table and histogram

%deputy_rtn 열어서 각 1:60:end마다 100개 반복에 대해 점수 min, max, med, avg, std 구하기 
T = table(t', mins, meds, means, maxs, ...
    'VariableNames',{'t_sec','min','median','mean','max'});
figure; hold on; grid on; box on
plot(t,mins,'-','DisplayName','Min','Color','blue');
%plot(t,meds,'-','DisplayName','median')
plot(t,means,'-','DisplayName','Mean','Color','red')
%plot(t,maxs,'-','DisplayName','max')
xlabel('Time (s)'); ylabel('$d_{\min}$ (m)','Interpreter','latex');
xlim([100 300]);
set(gca, 'YScale', 'log');  ylim([0.01 50]);
legend show
legend('Location','best')


%% statistic
idx = 1:Niter;
t   = idx(:);  
mins = zeros(numel(idx),1);
meds = zeros(numel(idx),1);
means = zeros(numel(idx),1);
maxs = zeros(numel(idx),1);
stds  = zeros(numel(idx),1);  % 추가
N_oob = zeros(numel(idx),1);  % 추가: 범위 밖 개수
N_nan = zeros(numel(idx),1);  % 추가: 최종 NaN 개수

xmin=rangeSDF(1); xmax=rangeSDF(2);
ymin=rangeSDF(3); ymax=rangeSDF(4);
zmin=rangeSDF(5); zmax=rangeSDF(6);

xvec = xmin:h:xmax;   % 길이 Nx
yvec = ymin:h:ymax;   % 길이 Ny
zvec = zmin:h:zmax;   % 길이 Nz

F = griddedInterpolant({xvec, yvec, zvec}, SDF, 'linear', 'none'); 

useGrid = ~isempty(SDF);
if useGrid, [Nx,Ny,Nz] = size(SDF); end
for k = 1:numel(idx)
    timeIdx = 1:60:size(deputy_rtn,2);
    X = squeeze(deputy_rtn(idx(k),timeIdx,1)); X = X(:);
    Y = squeeze(deputy_rtn(idx(k),timeIdx,2)); Y = Y(:);
    Z = squeeze(deputy_rtn(idx(k),timeIdx,3)); Z = Z(:);
    Nt_k = numel(X);
   
    ix = 1 + round((X - xmin)/h);
    iy = 1 + round((Y - ymin)/h);
    iz = 1 + round((Z - zmin)/h);
    in = X>=xmin & X<=xmax & Y>=ymin & Y<=ymax & Z>=zmin & Z<=zmax;
   
    phi = nan(Nt_k,1);

    %in = in & ix>=1 & ix<=Nx & iy>=1 & iy<=Ny & iz>=1 & iz<=Nz;
    ii = find(in);
    if ~isempty(ii)
        %lin = sub2ind([Nx,Ny,Nz], ix(ii), iy(ii), iz(ii));  % 벡터화 대입
        %phi(ii) = SDF(lin);
        phi(ii) = F(X(ii), Y(ii), Z(ii));
    end


    % 집계
    v = phi(~isnan(phi));
    if isempty(v)
        mins(k)=NaN; meds(k)=NaN; means(k)=NaN; maxs(k)=NaN; stds(k)=NaN;
    else
        mins(k)=min(v);  meds(k)=median(v);  means(k)=mean(v);
        maxs(k)=max(v);  stds(k)=std(v,0);
    end
    N_oob(k) = sum(~in);          % 범위 밖(조회 불가로 간주)
    N_nan(k) = sum(isnan(phi));   % 최종 NaN 총계(= 범위 밖 + 평가실패)
end
%%
% T = table(t(:), mins, meds, means, maxs, stds, ...
%     'VariableNames', {'iter','min','median','mean','max','std'});

figure; hold on; grid on; box on

% 주요 통계선
plot(t, mins,  '-', 'DisplayName','Min',  'Color','blue');   % blue
plot(t, means, '-', 'DisplayName','Mean', 'Color','red');% red

% 1σ 밴드 추가 (mean ± std)
% fill([t; flipud(t)], [means-stds; flipud(means+stds)], ...
%      [0.9 0.7 0.7], 'FaceAlpha',0.3, 'EdgeColor','none', ...
%      'DisplayName','±1σ range');

xlabel('Iteration number');
ylabel('$d_{\min}$ (m)','Interpreter','latex');
set(gca, 'YScale','log');
ylim([0.01 50]);
xlim([min(t) max(t)]);
legend('Location','best');
mean_mins = mean(mins)
dev_mins  = std(mins)
med_mins  = med(mins)
