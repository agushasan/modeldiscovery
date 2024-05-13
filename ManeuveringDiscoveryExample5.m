%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 25;
dt  = 0.001;
t   = dt:dt:tf;

n = 3;
r = 36;

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
C = [0 0 0 1 0 0;0 0 0 0 1 0; 0 0 0 0 0 1];

%% state initialization
x        = [0;0;0;0;0;0];
y        = [0;0;0];
vhat     = [0;0;0];
thetahat = zeros(r,1);
 
%% known paramaters
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
m33 = M(3,3);
mt  = m22*m33-m23*m32;

alpha1 = dt*(m-Yvd)/(m-Xud);
alpha2 = dt*(m*xg-Yrd)/(m-Xud);
alpha3 = (-dt*(Iz-Nrd)*(m-Xud)/mt)+(dt*(m*xg-Yrd)*(m*xg-Yrd)/mt);
alpha4 = dt*(m*xg-Yrd)*(Xud-Yvd)/mt;
alpha5 = (dt*(m*xg-Nvd)*(m-Xud)/mt)-(dt*(m-Yvd)*(m*xg-Yrd)/mt);
alpha6 = -dt*(m-Yvd)*(Xud-Yvd)/mt;

beta1  = dt*Xu/(m-Xud);
beta2  = dt*Xuu/(m-Xud);
beta3  = (dt*(Iz-Nrd)*Yv/mt)-(dt*(m*xg-Yrd)*Nv/mt);
beta4  = (dt*(Iz-Nrd)*Yr/mt)-(dt*(m*xg-Yrd)*Nr/mt);
beta5  = dt*(Iz-Nrd)*Yvv/mt;
beta6  = -dt*(m*xg-Yrd)*Nrr/mt;
beta7  = (dt*(m-Yvd)*Nv/mt)-(dt*(m*xg-Nvd)*Yv/mt);
beta8  = (dt*(m-Yvd)*Nr/mt)-(dt*(m*xg-Nvd)*Yr/mt);
beta9  = -dt*(m*xg-Nvd)*Yvv/mt;
beta10 = dt*(m-Yvd)*Nrr/mt;

%% initial control inputs
%u     = [40 10 1]';
u     = [5 10 0]';

%% for plotting
uArray          = [];
xArray          = [];
yArray          = [];
vhatArray       = [];
thetahatArray   = [];

%% Initialization for estimator

Pplus       = 100000*eye(n);
QF          = 10*eye(n);
RF          = 10*eye(n);
a           = 0.999;
UpsilonPlus = 0*zeros(n,r);
S           = 1*eye(r);
lambda      = 0.9999;

%% simulation
for i=1:(tf/dt)
    
    %u     = [40*sin(i*dt) 10*cos(i*dt) 1*cos(i*dt)]';

    uArray         = [uArray u];
    xArray         = [xArray x];
    yArray         = [yArray y];
    vhatArray      = [vhatArray vhat];
    thetahatArray  = [thetahatArray thetahat]; 

    Cvv = [alpha1*x(5)*x(6)+alpha2*x(6)^2;alpha3*x(4)*x(6)+alpha4*x(4)*x(5);alpha5*x(4)*x(6)+alpha6*x(4)*x(5)];
    Dvv = [beta1*x(4)+beta2*abs(x(4))*x(4);beta3*x(5)+beta4*x(6)+beta5*abs(x(5))*x(5)+beta6*abs(x(6))*x(6);beta7*x(5)+beta8*x(6)+beta9*abs(x(5))*x(5)+beta10*abs(x(6))*x(6)];

    x = x+[dt*(cos(x(3))*x(4)-sin(x(3))*x(5));dt*(sin(x(3))*x(4)+cos(x(3))*x(5));dt*x(6);Cvv+Dvv]+B*u;
    y = C*x;

%    Phi = [vhat(2)*vhat(3) vhat(3)^2 vhat(1) abs(vhat(1))*vhat(1) 0 0 0 0 0 0 0 0 0 0 0 0;
%          0 0 0 0 vhat(1)*vhat(3) vhat(1)*vhat(2) vhat(2) vhat(3) abs(vhat(2))*vhat(2) abs(vhat(3))*vhat(3) 0 0 0 0 0 0;
%          0 0 0 0 0 0 0 0 0 0 vhat(1)*vhat(3) vhat(1)*vhat(2) vhat(2) vhat(3) abs(vhat(2))*vhat(2) abs(vhat(3))*vhat(3)];

%     Phi = [vhat(1) vhat(2) vhat(3) vhat(1)^2 vhat(2)^2 vhat(3)^2 vhat(1)^3 vhat(2)^3 vhat(3)^3 zeros(1,18);
%            zeros(1,9) vhat(1) vhat(2) vhat(3) vhat(1)^2 vhat(2)^2 vhat(3)^2 vhat(1)^3 vhat(2)^3 vhat(3)^3 zeros(1,9);
%            zeros(1,18) vhat(1) vhat(2) vhat(3) vhat(1)^2 vhat(2)^2 vhat(3)^2 vhat(1)^3 vhat(2)^3 vhat(3)^3];

    Phi = [vhat(1) vhat(2) vhat(3) vhat(1)*vhat(2) vhat(1)*vhat(3) vhat(2)*vhat(3) vhat(1)^2 vhat(2)^2 vhat(3)^2 vhat(1)^3 vhat(2)^3 vhat(3)^3 zeros(1,24);
           zeros(1,12) vhat(1) vhat(2) vhat(3) vhat(1)*vhat(2) vhat(1)*vhat(3) vhat(2)*vhat(3) vhat(1)^2 vhat(2)^2 vhat(3)^2 vhat(1)^3 vhat(2)^3 vhat(3)^3 zeros(1,12);
           zeros(1,24) vhat(1) vhat(2) vhat(3) vhat(1)*vhat(2) vhat(1)*vhat(3) vhat(2)*vhat(3) vhat(1)^2 vhat(2)^2 vhat(3)^2 vhat(1)^3 vhat(2)^3 vhat(3)^3];

    % Estimation using adaptive KF
    Pmin  = Pplus+QF;
    Sigma = Pmin+RF;
    KF    = Pmin*inv(Sigma);
    Pplus = (eye(n)-KF)*Pmin;
     
    ytilde = y-(vhat+dt*inv(M)*(u)+Phi*thetahat);
    QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
    RF    = a*RF + (1-a)*(ytilde*ytilde'+Pmin);
 
    Upsilon = (eye(n)-KF)*(eye(n))*UpsilonPlus+(eye(n)-KF)*Phi;
    Omega   = eye(n)*UpsilonPlus+Phi;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma1  = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    thetahat  = thetahat + Gamma1*ytilde;
    vhat      = vhat+dt*inv(M)*(u)+Phi*thetahat+KF*ytilde+Upsilon*Gamma1*ytilde;

end

Temp1 = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetahatArray(7,:);thetahatArray(13,:)];
Temp2 = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetahatArray(8,:);thetahatArray(14,:)];

figure(1)
plot(xArray(1,:),xArray(2,:), ':b', 'LineWidth', 6)
hold on;
plot(xArray(1,1),xArray(2,1), 'ok', 'LineWidth', 20)
hold on;
plot(xArray(1,end),xArray(2,end), 'or', 'LineWidth', 20)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
legend('trajectory','start','end','FontSize',48)
xlabel('x (m)','FontSize',48)
ylabel('y (m)','FontSize',48)

figure(2)
subplot(3,1,1)
plot(t,yArray(1,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vhatArray(1,:), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('u [m/s]','FontSize',48)
subplot(3,1,2)
plot(t,yArray(2,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vhatArray(2,:), ':g', 'LineWidth', 6)
grid on;
grid minor;
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
ylabel('v [m/s]','FontSize',48)
legend('measured','estimated','FontSize',48)
subplot(3,1,3)
plot(t,yArray(3,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vhatArray(3,:), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('r [rad/s]','FontSize',48)
xlabel('time (s)','FontSize',48)

figure(3)
subplot(6,1,1)
plot(t,thetahatArray(1,:)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('\theta_1','FontSize',48)
subplot(6,1,2)
plot(t,thetahatArray(2,:)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('\theta_2','FontSize',48)
subplot(6,1,3)
plot(t,thetahatArray(3,:)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('\theta_3','FontSize',48)
subplot(6,1,4)
plot(t,thetahatArray(4,:)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('\theta_4','FontSize',48)
subplot(6,1,5)
plot(t,thetahatArray(5,:)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('\theta_5','FontSize',48)
subplot(6,1,6)
plot(t,thetahatArray(6,:)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('\theta_6','FontSize',48)
xlabel('time (s)','FontSize',48)

% figure(4)
% subplot(3,1,1)
% plot(t,thetahatArray(1,:)/dt, ':g', 'LineWidth', 6)
% set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
% grid on;
% grid minor;
% ylabel('\theta_1','FontSize',48)
% subplot(3,1,2)
% plot(t,thetahatArray(2,:)/dt, ':g', 'LineWidth', 6)
% set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
% grid on;
% grid minor;
% ylabel('\theta_2','FontSize',48)
% subplot(3,1,3)
% plot(t,thetahatArray(3,:)/dt, ':g', 'LineWidth', 6)
% set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
% grid on;
% grid minor;
% ylabel('\theta_3','FontSize',48)
% xlabel('time (s)','FontSize',48)

thetahatArray(:,end)/dt