%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 10;
dt  = 0.001;
t   = dt:dt:tf;

n = 3;
r = 16;

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
x        = [3;3;0;0;0;0];
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

Pplus       = 10000000000*eye(n);
QF          = 0.0001*eye(n);
RF          = 0.0001*eye(n);
a           = 0.999;
UpsilonPlus = 0*zeros(n,r);
S           = 10000000000*eye(r);
lambda      = 0.999;

%% simulation
for i=1:(tf/dt)
    
    u     = [40*sin(i*dt) 10*cos(i*dt) 1*cos(i*dt)]';

    uArray         = [uArray u];
    xArray         = [xArray x];
    yArray         = [yArray y];
    vhatArray      = [vhatArray vhat];
    thetahatArray  = [thetahatArray thetahat]; 

    Cvv = [alpha1*x(5)*x(6)+alpha2*x(6)^2;alpha3*x(4)*x(6)+alpha4*x(4)*x(5);alpha5*x(4)*x(6)+alpha6*x(4)*x(5)];
    Dvv = [beta1*x(4)+beta2*abs(x(4))*x(4);beta3*x(5)+beta4*x(6)+beta5*abs(x(5))*x(5)+beta6*abs(x(6))*x(6);beta7*x(5)+beta8*x(6)+beta9*abs(x(5))*x(5)+beta10*abs(x(6))*x(6)];

    x = x+[dt*(cos(x(3))*x(4)-sin(x(3))*x(5));dt*(sin(x(3))*x(4)+cos(x(3))*x(5));dt*x(6);Cvv+Dvv]+B*u;
    y = C*x+0.01*randn(3,1);

    Phi = [vhat(2)*vhat(3) vhat(3)^2 vhat(1) abs(vhat(1))*vhat(1) 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 vhat(1)*vhat(3) vhat(1)*vhat(2) vhat(2) vhat(3) abs(vhat(2))*vhat(2) abs(vhat(3))*vhat(3) 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 vhat(1)*vhat(3) vhat(1)*vhat(2) vhat(2) vhat(3) abs(vhat(2))*vhat(2) abs(vhat(3))*vhat(3)];

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
subplot(3,2,1)
plot(t,yArray(1,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vhatArray(1,:), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('u [m/s]','FontSize',48)
subplot(3,2,2)
plot(t,yArray(1,:), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vhatArray(1,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0.14 0.44])
xlim([5.5 6.5])
subplot(3,2,3)
plot(t,yArray(2,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vhatArray(2,:), ':g', 'LineWidth', 6)
grid on;
grid minor;
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
ylabel('v [m/s]','FontSize',48)
legend('measured','estimated','FontSize',48)
subplot(3,2,4)
plot(t,yArray(2,:), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vhatArray(2,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0.14 0.25])
xlim([0.53 1.5])
subplot(3,2,5)
plot(t,yArray(3,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vhatArray(3,:), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('r [rad/s]','FontSize',48)
xlabel('time (s)','FontSize',48)
subplot(3,2,6)
plot(t,yArray(3,:), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vhatArray(3,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 1])
ylim([0 0.14])
xlabel('time (s)','FontSize',48)

figure(2)
subplot(4,2,1)
plot(t,Xu*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),m11*thetahatArray(3,1:100:end)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Xu-4 Xu+4]);
ylabel('X_u','FontSize',48)
subplot(4,2,2)
plot(t,Xuu*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),m11*thetahatArray(4,1:100:end)/dt, ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('X_{uu}','FontSize',48)
legend('true parameter','estimated parameter','FontSize',36)
ylim([Xuu-4 Xuu+4]);
subplot(4,2,3)
plot(t,Nv*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp1(1,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Nv-4 Nv+4]);
ylabel('N_v','FontSize',48)
subplot(4,2,4)
plot(t,Yv*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp1(2,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('Y_v','FontSize',48)
ylim([Yv-4 Yv+4]);
subplot(4,2,5)
plot(t,Nr*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp2(1,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Nr-4 Nr+4]);
ylabel('N_r','FontSize',48)
subplot(4,2,6)
plot(t,Yr*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp2(2,1:100:end), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('Y_r','FontSize',48)
ylim([Yr-4 Yr+4]);
subplot(4,2,7)
plot(t,Yvv*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),mt*thetahatArray(9,1:100:end)/(dt*m33), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlabel('time (s)','FontSize',48)
ylim([Yvv-4 Yvv+4]);
ylabel('Y_{vv}','FontSize',48)
subplot(4,2,8)
plot(t,Nrr*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t(1:100:end),-mt*thetahatArray(10,1:100:end)/(dt*m23), ':g', 'LineWidth', 6)
set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('N_{rr}','FontSize',48)
xlabel('time (s)','FontSize',48)
ylim([Nrr-4 Nrr+4]);
