%2024 08 20

PP=generateRandomPointsOnSurface(403)+randn(403,3)*0.034;

PP2=PP(PP(:,3)>0.5,:);
PP3=PP(PP(:,1)>0.80|PP(:,1)<-.70,:);
PP4=PP(PP(:,2)>0.72&PP(:,3)>0,:);
PP5=PP(PP(:,2)<-0.73&PP(:,3)<0,:);
PP6=PP(PP(:,2)>0.652&PP(:,3)<0&PP(:,1)>0,:);

q= [cos(10/180*pi) 0 sin(10/180*pi)*1/sqrt(5) sin(10/180*pi)*2/sqrt(5)];
rotationM = quat2rotm(q);
R_PP2 = PP * rotationM';
 
%{
figure();
hold on
scatter3(PP2(:,1),PP2(:,2),PP2(:,3),'red');
scatter3(PP3(:,1),PP3(:,2),PP3(:,3),'b');
scatter3(PP4(:,1),PP4(:,2),PP4(:,3),'green');
scatter3(PP5(:,1),PP5(:,2),PP5(:,3),'black');
scatter3(PP6(:,1),PP6(:,2),PP6(:,3),'y');
hold off
%}
Ans1=regressionFourthOrder([PP2;PP4])
Ans2=regressionFourthOrder([PP2])
Ans3=regressionFourthOrder([PP])
Ans4=regressionFourthOrder([R_PP2])
f=@(x,y,z) Ans2(1)*x.^4+  Ans2(2)*y.^4+    Ans2(3)*z.^4+    Ans2(4)*x.^2*y.^2   +Ans2(5)*x.^2*z.^2   +Ans2(6)*y.^2*z.^2 ...
   +Ans2(7)*(x.^3).*y    +Ans2(8)*(x.^3).*z    +Ans2(9)*(y.^3).*x    +Ans2(10)*(y.^3).*z    -Ans2(11)*(z.^3).*x    +Ans2(12)*(z.^3).*y -1 ;
g=@(x,y,z) Ans1(1)*x.^4+  Ans1(2)*y.^4+    Ans1(3)*z.^4+    Ans1(4)*x.^2*y.^2   +Ans1(5)*x.^2*z.^2   +Ans1(6)*y.^2*z.^2 ...
   +Ans1(7)*(x.^3).*y    +Ans1(8)*(x.^3).*z    +Ans1(9)*(y.^3).*x    +Ans1(10)*(y.^3).*z    -Ans1(11)*(z.^3).*x    +Ans1(12)*(z.^3).*y -1 ;
h=@(x,y,z) Ans3(1)*x.^4+  Ans3(2)*y.^4+    Ans3(3)*z.^4+    Ans3(4)*x.^2*y.^2   +Ans3(5)*x.^2*z.^2   +Ans3(6)*y.^2*z.^2 ...
   +Ans3(7)*(x.^3).*y    +Ans3(8)*(x.^3).*z    +Ans3(9)*(y.^3).*x    +Ans3(10)*(y.^3).*z    +Ans3(11)*(z.^3).*x    +Ans3(12)*(z.^3).*y -1 ;
hr=@(x,y,z) Ans4(1)*x.^4+  Ans4(2)*y.^4+    Ans4(3)*z.^4+    Ans4(4)*x.^2*y.^2   +Ans4(5)*x.^2*z.^2   +Ans4(6)*y.^2*z.^2 ...
   +Ans4(7)*(x.^3).*y    +Ans4(8)*(x.^3).*z    +Ans4(9)*(y.^3).*x    +Ans4(10)*(y.^3).*z    +Ans3(11)*(z.^3).*x    +Ans4(12)*(z.^3).*y -1 ;


figure;
hold on
fimplicit3(g,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
scatter3(PP2(:,1),PP2(:,2),PP2(:,3),'b');
scatter3(PP4(:,1),PP4(:,2),PP4(:,3),'g');
hold off

figure;
hold on
fimplicit3(f,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
scatter3(PP2(:,1),PP2(:,2),PP2(:,3),'b');

hold off




figure;
hold on
%fimplicit3(h,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
fimplicit3(hr,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
scatter3(PP(:,1),PP(:,2),PP(:,3),'g');
scatter3(R_PP2(:,1),R_PP2(:,2),R_PP2(:,3),'b');

hold off