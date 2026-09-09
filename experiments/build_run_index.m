%% build_run_index.m — test_ridge/ 런 폴더 인덱스 생성
%
% 각 런이 어떤 축으로 돌았는지 test_ridge/RUNS.csv 로 모은다.
%   run_config.txt 가 있으면 그것을 읽고, 없으면(구버전 런) 결과 CSV 의
%   컬럼 고유값에서 축을 역산한다.
%
% 실행: run('experiments/build_run_index.m')

root   = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,'test_ridge');
d      = dir(outDir);
d      = d([d.isdir] & ~ismember({d.name},{'.','..'}));

R = table();
for i = 1:numel(d)
    f = fullfile(outDir, d(i).name);
    row = struct('folder',d(i).name,'kind',"",'stamp',"",'fold',"",'mode',"", ...
        'nRows',0,'nCells',0,'nSeed',0,'order',"",'tau',"",'lambda',"", ...
        'shape',"",'eps',"",'occ',"",'maxIter',"",'otype',"",'oracle',"",'note',"");

    % 종류: 폴더 이름 접두사
    if     startsWith(d(i).name,'sec32_single'),  row.kind = "single";
    elseif startsWith(d(i).name,'sec32_twoobj'),  row.kind = "twoobj";
    elseif startsWith(d(i).name,'sec32_diag'),    row.kind = "diag";
    elseif startsWith(d(i).name,'sec32_compare'), row.kind = "compare";
    else,                                          row.kind = "other";
    end
    tok = regexp(d(i).name,'(\d{8}_\d{4,6})','tokens','once');
    if ~isempty(tok), row.stamp = string(tok{1}); end
    if contains(d(i).name,'_fold'), row.fold = "eval(101-120)"; else, row.fold = "tuning(1-20)"; end
    if     contains(d(i).name,'_quick'), row.mode = "quick";
    elseif contains(d(i).name,'_dev'),   row.mode = "dev";
    else,                                row.mode = "";
    end

    csv = fullfile(f,'ransac_outlier_results.csv');
    if ~isfile(csv), csv = fullfile(f,'ransac_twoobj_results.csv'); end
    if isfile(csv)
        T = readtable(csv);
        row.nRows = height(T);
        row.nSeed = numel(unique(T.seed));
        keys = intersect({'shape','eps','lambda','maxIter','otype','occ','order','tau'}, ...
                         T.Properties.VariableNames);
        row.nCells = height(unique(T(:,keys)));
        row.order   = lv(T,'order');   row.tau     = lv(T,'tau');
        row.lambda  = lv(T,'lambda');  row.shape   = lv(T,'shape');
        row.eps     = lv(T,'eps');     row.occ     = lv(T,'occ');
        row.maxIter = lv(T,'maxIter'); row.otype   = lv(T,'otype');
        if ismember('cd_all_oracle',T.Properties.VariableNames)
            row.oracle = string(tern(any(~isnan(T.cd_all_oracle)),'ON','OFF'));
        end
    else
        row.note = "CSV 없음 (그림 전용 또는 중단)";
    end

    % run_config.txt 가 있으면 mode 를 그쪽 값으로 덮어쓴다 (신뢰도 높음)
    cfg = fullfile(f,'run_config.txt');
    if isfile(cfg)
        txt = string(fileread(cfg));
        m = regexp(txt,'mode = (\w+)','tokens','once');
        if ~isempty(m), row.mode = string(m{1}); end
    end

    R = [R; struct2table(row,'AsArray',true)]; %#ok<AGROW>
end

R = sortrows(R,'stamp','descend');
writetable(R, fullfile(outDir,'RUNS.csv'));
fprintf('%d개 런 → %s\n', height(R), fullfile(outDir,'RUNS.csv'));
disp(R(1:min(12,height(R)), {'folder','kind','mode','order','tau','lambda','occ','maxIter','nCells'}));

%% ---- local ----
function s = lv(T, name)
% 컬럼의 고유값을 짧은 문자열로. 값이 많으면 범위로 줄인다.
    if ~ismember(name, T.Properties.VariableNames), s = ""; return; end
    v = T.(name);
    if iscell(v) || isstring(v)
        s = strjoin(string(unique(v))', "|"); return;
    end
    u = unique(v);
    if numel(u) > 6
        s = sprintf('%g..%g (%d)', min(u), max(u), numel(u));
    else
        s = strjoin(arrayfun(@(x) string(num2str(x,'%g')), u)', "|");
    end
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
