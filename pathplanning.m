n  = sqrt(3.986*10^5/(6380+650)^3); 

D1_new = 0;
C2_new = -30;
C3_new = 0;
C4_new = -0;

x0_new = D1_new + C2_new;
y0_new = 2*C3_new+C4_new;
vx0_new = C3_new*n;
vy0_new = -3/2*n*D1_new-2*C2_new;

%%

% (2) 시간 범위 설정 (0 ~ 10초를 100개 구간으로 나눔)
t = linspace(0, 5400, 50);

% (3) x(t), y(t) 계산
x = D1_new + C2_new*cos(n*t) + C3_new*sin(n*t);
y = -3/2*n*D1_new.*t - 2*C2_new*sin(n*t) + 2*C3_new*cos(n*t) + C4_new;



% (5) 매개변수(Parametric) 플롯
figure;
plot(y, x, 'LineWidth', 1.2);
set(gca,'XDir','reverse')
xlabel('Y (m)','Fontsize',16);
ylabel('X (m)','Fontsize',16);

xlim([-200 200]);
ylim([-50 50]);
grid on;

disp([D1_new,C2_new,C3_new,C4_new]);
save('Path_new.mat',"D1_new","C2_new","C3_new","C4_new"); 
%%
x_max = -20;
x_min = -40;

D1_old = (x_max+x_min)/2;
C_norm = (x_max-x_min)/2;

x0_old = x0_new;
y0_old = y0_new;

C2_old = x0_old - D1_old;
C3_old = sqrt(C_norm^2-C2_old^2);
C4_old = y0_old - 2*C3_old;

t = linspace(-4000, 0, 100);
% (3) x(t), y(t) 계산
x = D1_old + C2_old*cos(n*t) + C3_old*sin(n*t);
y = -3/2*n*D1_old.*t - 2*C2_old*sin(n*t) + 2*C3_old*cos(n*t) + C4_old;



% (5) 매개변수(Parametric) 플롯
hold on;
plot(y, x, 'LineWidth', 1.2);
set(gca,'XDir','reverse')

hold off;
axis equal;
legend('Hovering Trajectory','Approach Trajectory')

disp([D1_old,C2_old,C3_old,C4_old]);

save('Path_old.mat',"D1_old","C2_old","C3_old","C4_old"); 