
clear; clc;
lambda=3/17.5;              %����
rho=3*lambda;           %���㵽ԭ�����
phi0=50*2*pi/360;         %�����
k=2*pi/lambda;           %����
phi=0:0.002:2*pi;           %����
n=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=phi-phi0;
beta2=phi+phi0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�缫��ƽ�沨ɢ�䳡:
Ud=f(beta1,n,rho,k)+f(-beta1,n,rho,k)-(f(beta2,n,rho,k)+f(-beta2,n,rho,k));
%�ż���ƽ�沨ɢ�䳡:
% Ud=f(beta1,n,rho,k)+f(-beta1,n,rho,k)+(f(beta2,n,rho,k)+f(-beta2,n,rho,k));
for jj=1:length(phi)
    Uo_i(jj)=exp(j*k*rho*cos(beta1(jj)));
    Uo_r(jj)=exp(j*k*rho*cos(beta2(jj)));
end
%�缫��ƽ�沨���ι�ѧ��:
for jj=1:length(phi)
Uo(jj)=Uo_i(jj)-Uo_r(jj);
end
%�ż���ƽ�沨���ι�ѧ��:
% for jj=1:length(phi)
% Uo(jj)=Uo_i(jj)+Uo_r(jj);
% end

plot(180/pi*phi,abs(Ud),'r--'),axis([0,360,0,2.5])
xlabel('�գ��㣩'),ylabel('��ǿ')
hold on
plot(180/pi*phi,abs(Uo+Ud),'b')
axis([0,360,0,2.5]),set(gca,'Xtick',[0:30:360],'Ytick',[0:0.5:2.5])
legend('���䳡','�ܳ�')
hold off
