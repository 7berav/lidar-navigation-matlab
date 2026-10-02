function [idx, Neff] = selectModes(U, nSkip, dim, minEff)
% 임베딩에 쓸 고유벡터 선택: 국소화된 모드를 건너뛴다
%
% 고유벡터 u (단위벡터) 의 유효 점 수(participation)
%       Neff = 1 / sum(u_i^4)
% 는 u 가 몇 개의 점에 퍼져 있는지를 나타낸다 (고르게 퍼지면 N, 한 점에 몰리면 1).
% Neff 가 작은 모드는 거의 끊어진 작은 조각에 국소화된 것으로, 크기 minEff 이상의
% 부품을 구분하는 데 쓸모가 없고 임베딩 차원만 차지한다.
%
% 입력
%   U      : N x kEig 고유벡터 (고유값 오름차순, embedSpectral 출력)
%   nSkip  : 앞에서 무조건 건너뛸 개수 (연결성분 수 = 영고유값 개수)
%   dim    : 뽑을 고유벡터 개수
%   minEff : 유효 점 수 하한 (기본 0 = 거르지 않음. 보통 minSize 와 같게)
%
% 출력
%   idx  : 선택된 열 번호 (dim 개. 모자라면 있는 만큼)
%   Neff : 1 x kEig 유효 점 수

if nargin < 4 || isempty(minEff), minEff = 0; end
Neff = 1 ./ sum(U.^4, 1);
cand = (nSkip+1):size(U,2);
cand = cand(Neff(cand) >= minEff);
idx  = cand(1:min(dim, numel(cand)));
end
