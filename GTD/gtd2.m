
clear; clc;
lambda=1/3;              %波长
rho=3*lambda;           %场点到原点距离
phi0=45*2*pi/360;         %入射角
k=2*pi/lambda;           %波数
phi=0:0.002:2*pi;           %步长
n=7/4;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=phi-phi0;
beta2=phi+phi0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%电极化平面波散射场:
% Ud=f(beta1,n,rho,k)+f(-beta1,n,rho,k)-(f(beta2,n,rho,k)+f(-beta2,n,rho,k));
%磁极化平面波散射场:
Ud=f(beta1,n,rho,k)+f(-beta1,n,rho,k)+(f(beta2,n,rho,k)+f(-beta2,n,rho,k));
for jj=1:length(phi)
if (phi(jj)>0)&(phi(jj)<pi-phi0)
    Uo_i(jj)=exp(j*k*rho*cos(beta1(jj)));
    Uo_r(jj)=exp(j*k*rho*cos(beta2(jj)));
elseif (phi(jj)>pi-phi0)&(phi(jj)<pi+phi0)
    Uo_i(jj)=exp(j*k*rho*cos(beta1(jj)));
    Uo_r(jj)=0;
else 
    Uo_i(jj)=0;
    Uo_r(jj)=0;
end
end
%电极化平面波几何光学场:
% for jj=1:length(phi)
% Uo(jj)=Uo_i(jj)-Uo_r(jj);
% end
%磁极化平面波几何光学场:
for jj=1:length(phi)
Uo(jj)=Uo_i(jj)+Uo_r(jj);
end
toc

plot(180/pi*phi,abs(Ud),'r--'),axis([0,360,0,2.5])
xlabel('φ（°）'),ylabel('场强')
hold on
plot(180/pi*phi,abs(Uo+Ud),'b')
axis([0,360,0,2.5]),set(gca,'Xtick',[0:30:360],'Ytick',[0:0.5:2.5])
legend('绕射场','总场')
hold off