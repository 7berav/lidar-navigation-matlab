function exportFuncToPython(FuncsC, outPyPath, argnames)
% FuncsC   : matlabFunction으로 만든 @(x,y,z) 형태의 함수핸들
% outPyPath: 생성할 .py 파일 경로 (예: 'poly_func.py')
% argnames : (선택) 인자명 셀배열. 기본 {'x','y','z'}

if nargin < 3
    argnames = {'x','y','z'};
end

% 1) 함수 본문 문자열화

s = func2str(FuncsC);


% 2) 헤더/줄바꿈/불필요 공백 처리
%    @(x,y,z) 제거, ... 줄바꿈 제거, 앞뒤 공백 제거
s = regexprep(s, '@\([^)]*\)', '');
s = regexprep(s, '\s*\.{3}\s*', '');
s = strtrim(s);

% 3) 연산자 치환 (순서 중요)
%    MATLAB: .^  .*  ./  ^   →  Python/NumPy: **  *  /  **
s = strrep(s, '.^', '**');
s = strrep(s, '.*', '*');
s = strrep(s, './', '/');
s = regexprep(s, '(?<!\*)\^(?!\*)', '**');   % 이미 바뀐 ** 제외하고 ^ → **

% 4) 함수명/상수 치환 (word boundary 사용)
pairs = { ...
  '\bsqrt\(',      'np.sqrt('; ...
  '\bsin\(',       'np.sin(';  '\bcos\(',     'np.cos(';  '\btan\(',   'np.tan('; ...
  '\basin\(',      'np.arcsin('; '\bacos\(',  'np.arccos('; '\batan2\(', 'np.arctan2('; '\batan\(', 'np.arctan('; ...
  '\bexp\(',       'np.exp(';   '\blog10\(',  'np.log10('; '\blog\(',   'np.log('; ...
  '\babs\(',       'np.abs(';   '\bsign\(',   'np.sign('; ...
  '\bmax\(',       'np.maximum('; '\bmin\(',  'np.minimum('; ...
  '\bpower\(',     'np.power(' };
for i=1:2:numel(pairs)
    s = regexprep(s, pairs{i}, pairs{i+1});
end

% 5) 상수/리터럴
s = regexprep(s, '\bpi\b', 'np.pi');
s = regexprep(s, '\bInf\b', 'np.inf');
s = regexprep(s, '\bNaN\b', 'np.nan');
s = regexprep(s, '\btrue\b', 'True');
s = regexprep(s, '\bfalse\b', 'False');

% 6) 안전 점검 (남은 MATLAB 전용 토큰이 있는지)
suspects = {'@',';',' end','^','.*','./','.\^'};
for k = 1:numel(suspects)
    if contains(s, suspects{k})
        warning('남은 MATLAB 토큰 감지: %s', suspects{k});
    end
end

% 7) Python 파일 생성
[fdir, fname, ~] = fileparts(outPyPath);
if ~isempty(fdir) && ~isfolder(fdir), mkdir(fdir); end

fid = fopen(outPyPath, 'w');
assert(fid>0, '파일을 열 수 없습니다: %s', outPyPath);

fprintf(fid, 'import numpy as np\n');
fprintf(fid, 'def %s(%s):\n', fname, strjoin(argnames, ', '));
fprintf(fid, '    return %s\n', s);
fclose(fid);

% 8) 식 자체(.txt)로도 저장(옵션, 디버깅용)
txtPath = fullfile(fdir, [fname '_expr.txt']);
fid = fopen(txtPath, 'w');
if fid>0
    fprintf(fid, '%s', s);
    fclose(fid);
end

fprintf('생성 완료: %s\n식 백업: %s\n', outPyPath, txtPath);
end