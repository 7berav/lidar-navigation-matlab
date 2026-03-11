syms x y z real

% 예제 다항식 (실제 수식 사용하세요)
f = x^6 + y^6 + z^6 + x^5*y + x^5*z + y^5*z + y^5*x + x^4*y^2 + x^2*y^4 + z^3*x^3;

[coeffs_raw, monomials_raw] = coeffs(expand(f), [x,y,z]);
num_terms = length(monomials_raw);

% 각 항의 전체 차수
total_degree = double(feval(symengine, 'degree', monomials_raw, [x,y,z]));

% 각 항에서의 개별 최대차수 (예: x^5*y는 max차수가 5)
degree_x = double(feval(symengine, 'degree', monomials_raw, x));
degree_y = double(feval(symengine, 'degree', monomials_raw, y));
degree_z = double(feval(symengine, 'degree', monomials_raw, z));

max_single_var_degree = max([degree_x; degree_y; degree_z]);

% 정렬 기준:
% 1순위: 전체차수 높은 순
% 2순위: 개별 최대 차수 높은 순
% 3순위: x차수 → y차수 → z차수 (옵션, 세부 동점처리)
sort_matrix = [total_degree; max_single_var_degree; degree_x; degree_y; degree_z]';

% 원하는 순서로 정렬 (내림차순)
[~, sortIdx] = sortrows(sort_matrix, [-1 -2 -3 -4 -5]);

% 재배열된 항과 계수
monomials_sorted = monomials_raw(sortIdx);
coeffs_sorted = coeffs_raw(sortIdx);

disp('사용자가 원하는 순서로 정렬된 항:');
disp(monomials_sorted.');

disp('대응하는 계수:');
disp(coeffs_sorted.');