Lx = 0.2; Ly = 0.3; Lz = 0.1;
x = Lx; y = Ly; z = Lz;

% 8개 꼭짓점 (좌표 중심 기준)
vertices = [
    -x -y -z;
     x -y -z;
     x  y -z;
    -x  y -z;
    -x -y  z;
     x -y  z;
     x  y  z;
    -x  y  z;
];
faces = [
    1 2 3 4;  % 아래
    5 6 7 8;  % 위
    1 2 6 5;  % 앞
    2 3 7 6;  % 오른쪽
    3 4 8 7;  % 뒤
    4 1 5 8;  % 왼쪽
];
figure;
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [0.95, 0.82, 0.5], 'FaceAlpha', 0.95, 'EdgeColor', 'k');
axis equal;
xlabel('X (m)',FontSize=16); ylabel('Y (m)',FontSize=16); zlabel('Z (m)',FontSize=16);
view(3);
xlim([-0.6 0.6]); ylim([-0.6 0.6]); zlim([-0.6 0.6]);
    set(gca, 'Color', 'w');
%%

filename = 'cubesat2.gif';
figure(9)
set(9, 'Color', 'w');
for i = 1:30
    clf
    rotm = quat2rotm(quater1(i,:));
    rotated_vertices = (rotm * vertices.').';  % 12×3
    offset = [0, 0, 0];
    moved_vertices = rotated_vertices - offset;

    patch('Vertices', moved_vertices, 'Faces', faces, ...
          'FaceColor', [0.95, 0.82, 0.5], 'EdgeColor', 'k', 'FaceAlpha', 0.95);
    set(gca, 'Color', 'w'); 

    axis equal; grid on;

    view(3); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    set(gca, 'Color', 'w');
    xlim([-0.6 0.6]); ylim([-0.6 0.6]); zlim([-0.6 0.6]);

    drawnow;

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end
%%
TermsB = homogeneTerm(order);
%fA = BetaA *  TermsB(:);
fB = Array_BetaB(131,2:29) *  TermsB(:);
%fA_shift =  subs(fA, [x, y, z], [x - DispA(1), y - DispA(2), z - DispA(3)])-1;
fB_shift =  subs(fB, [x, y, z], [x - Array_Disp(131,2), y - Array_Disp(131,3), z - Array_Disp(131,4)]) -1;


figure(4)
fimplicit3(fB_shift,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);

colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis [-1.2 1.2 -1.2 1.2 -1.2 1.2])
grid on
