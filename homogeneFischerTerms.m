function fisherTerms = homogeneFischerTerms(order)
% HOMOGENEFISCHERTERMS_LEGENDRE  고체 구면조화함수(Fischer 분해) 생성
%   입력: order (비음수 정수) — 원하는 최고차 ℓ
%   출력: fisherTerms — r^(order−ℓ)·H_ℓ^m(x,y,z) 를 세로 벡터로 반환
%             ℓ = order, order−2, …, 0  /  m = −ℓ…ℓ
%   ▷ 연관 르장드르 다항식은 legendreP(ℓ,m,x) 3-인수 버전이 MATLAB 내장에 없으므로
%     Rodrigues 공식으로 직접 생성한다.

% 1. 심볼릭 변수
    syms x y z t real                                % 직교 좌표
    r        = sqrt(x^2 + y^2 + z^2);              % 반경 r
    cosTheta = z / r;                              % cosθ = z/r
    phi      = atan2(y, x);                        % 방위각 φ

% 2. 결과 벡터 선할당
    levels     = order:-2:0;                       % ℓ = order, order−2, …
    numTerms   = sum(levels + 1);                % 총 항 개수 (각 ℓ마다 2ℓ+1)
    fisherTerms = sym(zeros(1,numTerms));          % 결과 초기화

    idx = 1;                                       % 결과 인덱스
    for l = levels
        for m = 0:l
            % 3. Rodrigues 공식으로 P_l^m(cosθ)
            %t  = cosTheta;                         % 치환용 변수 t
            Pl = 1/(2^l*factorial(l))*diff((t^2-1)^l, t, l);   % P_l(t)

            if m > 0
                Plm = (1 - t^2)^(m/2) * diff(Pl, t, m);        % P_l^m
            else
                Plm = Pl;
            end
            Plm_t = subs(Plm, t, cosTheta);
            % 4. 정규화 상수와 Y_l^m(θ,φ)
            K   = sqrt(factorial(l-abs(m))/factorial(l+abs(m))); % Condon–Shortley
            Ylm = K * Plm_t * exp(1i*m*phi);
            %Ylm = Plm_t * exp(1i*m*phi);
            % 5. 고체 구면조화 다항식 H_lm(x,y,z)
            Hxyz = simplify(expand(r^order * Ylm));
            % 6. 실수 기저: m=0 한 개, m>0 두 개 (cos, sin)
            if m == 0
                fisherTerms(idx) = simplify(real(Hxyz));
                idx = idx + 1;
            else
                %fisherTerms(idx)   = simplify(1*real(Hxyz));    % cos(mφ)
                fisherTerms(idx)   = simplify(sqrt(2)*real(Hxyz));    % cos(mφ)
                %fisherTerms(idx)   = fisherTerms(idx)/gcd(coeffs(fisherTerms(idx)));
                %fisherTerms(idx+1) = simplify(1*imag(Hxyz));    % sin(mφ)
                fisherTerms(idx+1) = simplify(sqrt(2)*imag(Hxyz));    % sin(mφ)
                %fisherTerms(idx+1) = fisherTerms(idx+1)/gcd(coeffs(fisherTerms(idx+1)));
                idx = idx + 2;
            end
        end
    end
end