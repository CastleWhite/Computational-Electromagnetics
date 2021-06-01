
clear; clc;
lambda=1/3;              %����
rho=3*lambda;           %���㵽ԭ�����
phi0=45*2*pi/360;         %�����
k=2*pi/lambda;           %����
phi=0:0.002:2*pi;           %����
n=7/4;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=phi-phi0;
beta2=phi+phi0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�缫��ƽ�沨ɢ�䳡:
% Ud=f(beta1,n,rho,k)+f(-beta1,n,rho,k)-(f(beta2,n,rho,k)+f(-beta2,n,rho,k));
%�ż���ƽ�沨ɢ�䳡:
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
%�缫��ƽ�沨���ι�ѧ��:
% for jj=1:length(phi)
% Uo(jj)=Uo_i(jj)-Uo_r(jj);
% end
%�ż���ƽ�沨���ι�ѧ��:
for jj=1:length(phi)
Uo(jj)=Uo_i(jj)+Uo_r(jj);
end
toc

plot(180/pi*phi,abs(Ud),'r--'),axis([0,360,0,2.5])
xlabel('�գ��㣩'),ylabel('��ǿ')
hold on
plot(180/pi*phi,abs(Uo+Ud),'b')
axis([0,360,0,2.5]),set(gca,'Xtick',[0:30:360],'Ytick',[0:0.5:2.5])
legend('���䳡','�ܳ�')
hold off