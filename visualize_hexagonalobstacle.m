H = 2;  % 높이
r = 1;  % 반지름 (단면의 크기)
theta = (0:5)*(2*pi/6);  % 0~300도, 6각형

% 아래 육각형 정점 (z = -H/2)
x1 = r * cos(theta);
y1 = r * sin(theta);
z1 = -H/2 * ones(1,6);

% 위 육각형 정점 (z = +H/2)
x2 = x1;
y2 = y1;
z2 = H/2 * ones(1,6);

% 총 12점 (각각 위아래 6점씩)
vertices = [x1', y1', z1'; x2', y2', z2'];  % (12×3)

faces = [
    1 2 8 7;  % 옆면 1
    2 3 9 8;  % 옆면 2
    3 4 10 9;
    4 5 11 10;
    5 6 12 11;
    6 1 7 12;  % 옆면 6
    1 2 3 4;
    4 5 6 1;% 아래면
    7 8 9 10;
    10 11 12 7];  % 윗면


figure;
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [0.95, 0.82, 0.5], 'FaceAlpha', 0.95, 'EdgeColor', 'k');
axis equal;
xlabel('X (m)',FontSize=16); ylabel('Y (m)',FontSize=16); zlabel('Z (m)',FontSize=16);
view(3);
axis ([-1.5 1.5 -1.5 1.5 -1.5 1.5]);

%%

filename = 'hexagonal_prism5.gif';
figure
for i = 261:297
    clf
    rotm = quat2rotm(quater0(i,:));
    rotated_vertices = (rotm * vertices.').';  % 12×3
    offset = [Orbit(i,2), Orbit(i,3), 0];
    moved_vertices = rotated_vertices - offset;

    patch('Vertices', moved_vertices, 'Faces', faces, ...
          'FaceColor', [0.95, 0.82, 0.5], 'EdgeColor', 'k', 'FaceAlpha', 0.95);
    set(gca, 'Color', 'w'); 

    axis equal; grid on;
    xlim([-25 25]);
    ylim([30 60]);
    zlim([-25 25]);
    view(3); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    set(gca, 'Color', 'w');
    drawnow;

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    if i == 261
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
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
