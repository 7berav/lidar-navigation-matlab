function points = generatePointsOnIcosahedron( N)
    phi = (1 + sqrt(5)) / 2;

    % 2) 정20면체의 12개 꼭짓점 좌표 정의
    V = [
        0,   1,   phi;
        0,  -1,   phi;
        0,   1,  -phi;
        0,  -1,  -phi;
        1,   phi,  0;
       -1,   phi,  0;
        1,  -phi,  0;
       -1,  -phi,  0;
        phi,  0,   1;
       -phi,  0,   1;
        phi,  0,  -1;
       -phi,  0,  -1
    ];

    F = [
        1, 9,  2;
        1, 2, 10;
        1, 10, 6;
        1, 6, 5;
        1, 5, 9;
    
        9, 5, 11;
        5, 6, 3;
        6, 10, 12;
        10, 2, 8;
        2, 9, 7;
    
        4, 11, 3;
        4, 3, 12;
        4, 12, 8;
        4, 8, 7;
        4, 7, 11;
    
        11, 7, 9;
        3, 11, 5;
        12, 3, 6;
        8, 12, 10;
        7, 8, 2
    ];

    for i = 1:size(V, 1)
        V(i, :) = V(i, :) / norm(V(i, :));
    end

    points = zeros(N, 3);
    num_faces = 20;

    % 삼각형 당 몇 개의 점을 만들지 결정
    points_per_face = ceil(N / num_faces);
    count = 0;

    for i = 1:num_faces
        v1 = V(F(i, 1), :);
        v2 = V(F(i, 2), :);
        v3 = V(F(i, 3), :);

        % 로컬 좌표계
        e1 = v2 - v1;
        e2 = v3 - v1;

        % 각 삼각형에서 점 찍기
        for j = 1:points_per_face
            % 삼각형 내부 균일분포 샘플링
            a = rand(); b = rand();
            if a + b > 1
                a = 1 - a;
                b = 1 - b;
            end
            % (a, b)는 배리센트릭 계수 → 3D 점 변환
            p = v1 + a * e1 + b * e2;

            count = count + 1;
            if count <= N
                points(count, :) = p;
            else
                return;
            end
        end
    end
end
