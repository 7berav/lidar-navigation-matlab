%% compare_outlier_types.m — 오염 유형 3종 정면 비교 (Sec 3.2)
%
% ransac_outlier_sweep.m(균일 / 배경평면) + ransac_two_object.m(구조적 2물체)의
% 최신 결과를 읽어 같은 ε 축에서 겹쳐 그린다.
%
% 오염 3종: 배경평면 / 균일 랜덤 / 2물체.
% F1 외에 fp_lev(오분류 점의 GT 표면까지 거리)와 손실 분해(cd_all - cd_all_oracle)를
% 같이 그린다.

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))
outRoot=fullfile(fileparts(mfilename('fullpath')),'..','..','test_ridge');
% quickMode=true 이면 _quick 산출물도 후보에 넣는다(스모크 테스트용).
if ~exist('quickMode','var') || isempty(quickMode), quickMode = false; end

Ts=readtable(fullfile(latestDir(outRoot,'sec32_single_',quickMode),'ransac_outlier_results.csv'));
Tt=readtable(fullfile(latestDir(outRoot,'sec32_twoobj_',quickMode),'ransac_twoobj_results.csv'));

baseLam = mode(Ts.lambda(Ts.lambda>0));
mbMaxIter = mode(Ts.maxIter);
selS = @(ot) strcmp(Ts.shape,'hexagon') & abs(Ts.lambda-baseLam)<1e-12 & ...
             Ts.maxIter==mbMaxIter & strcmp(Ts.otype,ot);
selT = strcmp(Tt.pair,'hexagon+cube') & strcmp(Tt.offset,'adjacent') & strcmp(Tt.scheme,'perobj');

% 실패 시드를 뺀 중앙값 (failRate 는 따로 본다)
medS = @(ot,e,v) median(Ts.(v)(selS(ot) & abs(Ts.eps-e)<1e-9 & Ts.fail==0),'omitnan');
medT = @(e,v)    median(Tt.(v)(selT     & abs(Tt.eps-e)<1e-9 & Tt.fail==0),'omitnan');
epsP = unique(Ts.eps(selS('plane')))';
epsU = unique(Ts.eps(selS('uniform')))';
epsT = unique(Tt.eps(selT))';

series = { 'background plane', epsP, @(e,v) medS('plane',e,v) ; ...
           'uniform',          epsU, @(e,v) medS('uniform',e,v) ; ...
           'structured (2obj)',epsT, @(e,v) medT(e,v) };

figDir=fullfile(outRoot,sprintf('sec32_compare_%s',datestr(now,'yyyymmdd_HHMM'))); mkdir(figDir);

fig=figure('Visible','off','Position',[50 50 1500 800]);
panels = { 'F1',      'F1^{inlier} (median)',      '오염 유형별 분류 성능' ; ...
           'precision','precision',                'precision 이 먼저 무너진다' ; ...
           'cd_all',  'CD\_all',                   '기하 오차' ; ...
           'fp_lev',  'fp\_lev (FP→GT 표면 거리)', '잘못 넣은 점이 표면에서 떨어진 거리' };
for pj=1:4
    subplot(2,3,pj); hold on; grid on;
    for si=1:size(series,1)
        ee=series{si,2}; f=series{si,3};
        y=arrayfun(@(e) f(e,panels{pj,1}),ee);
        if all(isnan(y)), continue; end
        plot(100*ee,y,'-o','LineWidth',2,'DisplayName',series{si,1});
    end
    xlabel('\epsilon [%]'); ylabel(panels{pj,2}); title(panels{pj,3});
    legend('Location','best','FontSize',7);
    if any(strcmp(panels{pj,1},{'F1','precision'})), ylim([0 1]); end
end

% 손실 분해: 분류 실패에 귀속되는 몫 (L3 − L2)
subplot(2,3,5); hold on; grid on;
for si=1:size(series,1)
    ee=series{si,2}; f=series{si,3};
    y=arrayfun(@(e) f(e,'cd_all') - f(e,'cd_all_oracle'),ee);
    if all(isnan(y)), continue; end
    plot(100*ee,y,'-o','LineWidth',2,'DisplayName',series{si,1});
end
yline(0,'k--'); xlabel('\epsilon [%]'); ylabel('CD\_all − CD\_all\_oracle');
title('분류 실패에 귀속되는 손실 (L3−L2)'); legend('Location','best','FontSize',7);

% 오라클 상한은 오염 유형과 무관해야 한다(정합성 확인용)
subplot(2,3,6); hold on; grid on;
for si=1:size(series,1)
    ee=series{si,2}; f=series{si,3};
    y=arrayfun(@(e) f(e,'cd_all_oracle'),ee);
    if all(isnan(y)), continue; end
    plot(100*ee,y,'-o','LineWidth',2,'DisplayName',series{si,1});
end
xlabel('\epsilon [%]'); ylabel('CD\_all\_oracle');
title('L2 상한 (오염 유형에 불변이어야 정상)'); legend('Location','best','FontSize',7);

sgtitle(sprintf('아웃라이어 유형 비교 (hexagon, \\lambda=%.0e)',baseLam));
exportgraphics(fig,fullfile(figDir,'fig4_outlier_type_comparison.png'),'Resolution',150);
savefig(fig,fullfile(figDir,'fig4_outlier_type_comparison.fig'));

%% ---- 콘솔 요약 ----
fprintf('\n=== 오염 유형 비교 (hexagon, λ=%.0e, 실패 제외 중앙값) ===\n', baseLam);
fprintf('%-5s | %-22s | %-22s | %-22s\n','eps','background plane','uniform','structured (2obj)');
fprintf('%-5s | %-22s | %-22s | %-22s\n','', 'F1/prec/CD/fp_lev','F1/prec/CD/fp_lev','F1/prec/CD/fp_lev');
for e = unique([epsP epsU epsT])
    cells = cell(1,3);
    for si=1:3
        f=series{si,3};
        if ismember(e, series{si,2})
            cells{si}=sprintf('%.2f/%.2f/%.2f/%.2f', f(e,'F1'), f(e,'precision'), f(e,'cd_all'), f(e,'fp_lev'));
        else
            cells{si}='-';
        end
    end
    fprintf('%-5.2f | %-22s | %-22s | %-22s\n', e, cells{1}, cells{2}, cells{3});
end
fprintf('그림 → %s\n',figDir);

function d=latestDir(outRoot,prefix,allowQuick)
    dd=dir(outRoot); nm={dd([dd.isdir]).name};
    if allowQuick
        hits=sort(nm(startsWith(nm,prefix)));                    % _quick 포함, 최신 우선
    else
        hits=sort(nm(startsWith(nm,prefix)&~contains(nm,'_quick')));
    end
    assert(~isempty(hits),'no run dir: %s',prefix);
    d=fullfile(outRoot,hits{end});
end
