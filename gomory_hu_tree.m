function [T] = gomory_hu_tree(G)
    N = numnodes(G);
    D = make_digraph_from_ugraph(G);
    parent = ones(N,1); parent(1)=0;
    cap = zeros(N,1);
    for s = 2:N
        t = parent(s);

        % --- maxflow 호출 ---
        [fval,~,Sraw] = maxflow(D, s, t);

        % --- [중요] Sraw 형식 통일 (논리 마스크 Smask) ---
        if islogical(Sraw)
            % graph일 때: 논리 벡터
            Smask = Sraw(:);
            if numel(Smask) < N, Smask(end+1:N) = false; end
        else
            % digraph일 때: 정수 인덱스 리스트
            Smask = false(N,1);
            Sraw = Sraw(:);
            Sraw = Sraw(Sraw>=1 & Sraw<=N);   % 방어적 클리핑
            Smask(Sraw) = true;
        end

        cap(s) = fval;

        % (옵션) t가 S쪽이면 GH 표준 절차에 따라 스와핑/재배치가 필요
        % 실용상 최소 구현: parent 재배치만 수행
        for v = 1:N
            if v ~= s && parent(v) == t && Smask(v)
                parent(v) = s;
            end
        end
        
        % 혹시라도 parent(s)==s가 되면 원래 t로 복원
        if parent(s) == s
            parent(s) = t;
        end
    end

    sT = (2:N)'; 
    tT = parent(2:N);
    wT = cap(2:N);
    T  = graph(sT, tT, wT, N);
end