%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 25;
dt  = 0.001;
t   = dt:dt:tf;

%% True parameters

m    = 23.8;
Iz   = 1.76;
xg   = 0.046;

Xud  = -2;
Yvd  = -10;
Yrd  = 0;
Nvd  = 0;
Nrd  = -1;

Xu   = -0.7225;
Xuu  = -1.3274;
Yv   = -0.8612;
Yvv  = -36.2823;
Yr   = 0.1079;
Nv   = 0.1052;
Nr  = -0.5;
Nrr = -1;

%% system description
M = [m-Xud 0 0;0 m-Yvd m*xg-Yrd;0 m*xg-Nvd Iz-Nrd];
B = dt*[zeros(3);inv(M)];

%% state initialization
x        = [0;0;0;0;0;0];
xD1      = [0;0;0;0;0;0];
xD2      = [0;0;0;0;0;0];
 
%% known paramaters
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
m33 = M(3,3);
mt  = m22*m33-m23*m32;

alpha1 = dt*(m-Yvd)/(m-Xud)
alpha2 = dt*(m*xg-Yrd)/(m-Xud)
alpha3 = (-dt*(Iz-Nrd)*(m-Xud)/mt)+(dt*(m*xg-Yrd)*(m*xg-Yrd)/mt)
alpha4 = dt*(m*xg-Yrd)*(Xud-Yvd)/mt
alpha5 = (dt*(m*xg-Nvd)*(m-Xud)/mt)-(dt*(m-Yvd)*(m*xg-Yrd)/mt)
alpha6 = -dt*(m-Yvd)*(Xud-Yvd)/mt

beta1  = dt*Xu/(m-Xud)
beta2  = dt*Xuu/(m-Xud)
beta3  = (dt*(Iz-Nrd)*Yv/mt)-(dt*(m*xg-Yrd)*Nv/mt)
beta4  = (dt*(Iz-Nrd)*Yr/mt)-(dt*(m*xg-Yrd)*Nr/mt)
beta5  = dt*(Iz-Nrd)*Yvv/mt
beta6  = -dt*(m*xg-Yrd)*Nrr/mt
beta7  = (dt*(m-Yvd)*Nv/mt)-(dt*(m*xg-Nvd)*Yv/mt)
beta8  = (dt*(m-Yvd)*Nr/mt)-(dt*(m*xg-Nvd)*Yr/mt)
beta9  = -dt*(m*xg-Nvd)*Yvv/mt
beta10 = dt*(m-Yvd)*Nrr/mt

%% initial control inputs
u     = [5 10 0]';

%% for plotting
uArray          = [];
xArray          = [];
xD1Array        = [];
xD2Array        = [];

%% Theta
Theta1 = [-0.0232 0.2468 0.5817 -0.0636 -0.3105 -0.1717 0.0208 -0.4074 -0.0823;
          0.2322 -0.2750 -0.0259 0.1514 -0.5004 -0.0582 0.0991 -0.2701 -0.1160;
          -1.4407 0.5636 -0.2305 -0.2972 -0.4264 0.3680 0.0248 -0.3594 -0.0810]';

Theta2 = [0.0245 0.1531 0.3997 -0.0797 0.0405 0.4023 -0.0401 -0.1952 -0.1158 0.0158 -0.2571 -0.0686;
          0.0747 -0.2356 -0.0552 -0.0565 -0.7483 0.1365 0.0239 -0.4826 -0.0209 0.0175 -0.3769 -0.0219;
          -0.9431 0.3997 -0.2291 -0.9278 0.1410 -0.0097 -0.1719 -0.2144 0.3278 0.0331 -0.1499 -0.1146]';

%% simulation
for i=1:(tf/dt)

    uArray         = [uArray u];
    xArray         = [xArray x];
    xD1Array       = [xD1Array xD1];
    xD2Array       = [xD2Array xD2];    

    Cvv = [alpha1*x(5)*x(6)+alpha2*x(6)^2;alpha3*x(4)*x(6)+alpha4*x(4)*x(5);alpha5*x(4)*x(6)+alpha6*x(4)*x(5)];
    Dvv = [beta1*x(4)+beta2*abs(x(4))*x(4);beta3*x(5)+beta4*x(6)+beta5*abs(x(5))*x(5)+beta6*abs(x(6))*x(6);beta7*x(5)+beta8*x(6)+beta9*abs(x(5))*x(5)+beta10*abs(x(6))*x(6)];

    x = x+[dt*(cos(x(3))*x(4)-sin(x(3))*x(5));dt*(sin(x(3))*x(4)+cos(x(3))*x(5));dt*x(6);Cvv+Dvv]+B*u;

    Phi1 = [xD1(4) xD1(5) xD1(6) xD1(4)^2 xD1(5)^2 xD1(6)^2 xD1(4)^3 xD1(5)^3 xD1(6)^3 zeros(1,18);
           zeros(1,9) xD1(4) xD1(5) xD1(6) xD1(4)^2 xD1(5)^2 xD1(6)^2 xD1(4)^3 xD1(5)^3 xD1(6)^3 zeros(1,9);
           zeros(1,18) xD1(4) xD1(5) xD1(6) xD1(4)^2 xD1(5)^2 xD1(6)^2 xD1(4)^3 xD1(5)^3 xD1(6)^3];

    Phi2 = [xD2(4) xD2(5) xD2(6) xD2(4)*xD2(5) xD2(4)*xD2(6) xD2(5)*xD2(6) xD2(4)^2 xD2(5)^2 xD2(6)^2 xD2(4)^3 xD2(5)^3 xD2(6)^3 zeros(1,24);
           zeros(1,12) xD2(4) xD2(5) xD2(6) xD2(4)*xD2(5) xD2(4)*xD2(6) xD2(5)*xD2(6) xD2(4)^2 xD2(5)^2 xD2(6)^2 xD2(4)^3 xD2(5)^3 xD2(6)^3 zeros(1,12);
           zeros(1,24) xD2(4) xD2(5) xD2(6) xD2(4)*xD2(5) xD2(4)*xD2(6) xD2(5)*xD2(6) xD2(4)^2 xD2(5)^2 xD2(6)^2 xD2(4)^3 xD2(5)^3 xD2(6)^3];

    xD1 = xD1+[dt*(cos(xD1(3))*xD1(4)-sin(xD1(3))*xD1(5));dt*(sin(xD1(3))*xD1(4)+cos(xD1(3))*xD1(5));dt*xD1(6);Phi1*Theta1*dt]+B*u;

    xD2 = xD2+[dt*(cos(xD2(3))*xD2(4)-sin(xD2(3))*xD2(5));dt*(sin(xD2(3))*xD2(4)+cos(xD2(3))*xD2(5));dt*xD2(6);Phi2*Theta2*dt]+B*u;

end

figure(1)
plot(xArray(1,:),xArray(2,:), ':b', 'LineWidth', 6)
hold on;
plot(xArray(1,1),xArray(2,1), 'ok', 'LineWidth', 20)
hold on;
plot(xArray(1,end),xArray(2,end), 'or', 'LineWidth', 20)
hold on;
plot(xD1Array(1,:),xD1Array(2,:), ':g', 'LineWidth', 6)
hold on;
plot(xD1Array(1,1),xD1Array(2,1), 'om', 'LineWidth', 20)
hold on;
plot(xD1Array(1,end),xD1Array(2,end), 'oc', 'LineWidth', 20)
hold on;
plot(xD2Array(1,:),xD2Array(2,:), ':y', 'LineWidth', 6)
hold on;
plot(xD2Array(1,1),xD2Array(2,1), 'om', 'LineWidth', 20)
hold on;
plot(xD2Array(1,end),xD2Array(2,end), 'oc', 'LineWidth', 20)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
legend('measurement','start','end','self-discovery 1','start','end','self-discovery 2','start','end','FontSize',48)
xlabel('x (m)','FontSize',48)
ylabel('y (m)','FontSize',48)

figure(2)
plot(xArray(1,:),xArray(2,:), ':b', 'LineWidth', 6)
hold on;
plot(xArray(1,1),xArray(2,1), 'ok', 'LineWidth', 20)
hold on;
plot(xArray(1,end),xArray(2,end), 'or', 'LineWidth', 20)
hold on;
plot(xD1Array(1,:),xD1Array(2,:), ':g', 'LineWidth', 6)
hold on;
plot(xD1Array(1,1),xD1Array(2,1), 'om', 'LineWidth', 20)
hold on;
plot(xD1Array(1,end),xD1Array(2,end), 'oc', 'LineWidth', 20)
hold on;
plot(xD2Array(1,:),xD2Array(2,:), ':y', 'LineWidth', 6)
hold on;
plot(xD2Array(1,1),xD2Array(2,1), 'om', 'LineWidth', 20)
hold on;
plot(xD2Array(1,end),xD2Array(2,end), 'oc', 'LineWidth', 20)
set(gca,'XTick',[], 'YTick', [],'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;

RMSExD1 = rmse(xArray(1,:),xD1Array(1,:))+rmse(xArray(2,:),xD1Array(2,:))
RMSExD2 = rmse(xArray(1,:),xD2Array(1,:))+rmse(xArray(2,:),xD2Array(2,:))
