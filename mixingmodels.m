%% 원래 가지고 있던 파일로 개형 만들기


E1=load("equation_total_margin.mat");
E2=load("equation_panel1_margin.mat");
E3=load("equation_panel2_margin.mat");
[coeffs1, terms1]=coeffs(E1.f3_translated_expanded, [x y z]);
[coeffs2, terms2]=coeffs(E2.f3_translated_expanded, [x y z]);
[coeffs3, terms3]=coeffs(E3.f3_translated_expanded, [x y z]);

filename = 'equations_output.txt';
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', equation2string(coeffs1, terms1));
fprintf(fid, '\n')
fprintf(fid, '%s\n', equation2string(coeffs2, terms2));
fprintf(fid, '\n')
fprintf(fid, '%s\n', equation2string(coeffs3, terms3));

fclose(fid);  % 파일 닫기

figure
hold on
fimplicit3(E1.f3_translated_expanded,[-3 3 -4 4 -2 2],'FaceColor', [0.99 0.75 0.12]);
fimplicit3(E2.f3_translated_expanded,[-2 2 -4 4 -2 2],'FaceColor', [0.05, 0.2, 0.5]);
fimplicit3(E3.f3_translated_expanded,[-2 2 -4 4 -2 2],'FaceColor', [0.05, 0.2, 0.5]);
hold off


colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
axis equal
grid on
view([1 1 1]);
xlim([-1 1]);  
ylim([-4 4]);
zlim([-2 2]);

figure
hold on
scatter3(PPm_total(:,1),PPm_total(:,2),PPm_total(:,3),3,PPm_total(:,3),'filled');
scatter3(PPm_panel1(:,1),PPm_panel1(:,2),PPm_panel1(:,3),1,PPm_panel1(:,3),'filled');
scatter3(PPm_panel2(:,1),PPm_panel2(:,2),PPm_panel2(:,3),1,PPm_panel2(:,3),'filled');
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
axis equal
grid on
view([1 1 1]);
xlim([-1 1]);  
ylim([-4 4]);
zlim([-2 2]);