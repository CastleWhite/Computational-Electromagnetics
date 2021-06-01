%计算 子函数
function re=f(beta,n,rho,k)
phi=0:0.002:2*pi;           %步长
N=200;
alpha=1+cos(beta);
for jj=1:length(phi)
t=sqrt(k*rho*alpha(jj)):0.01:N;
fun1=exp(-j*t.^2);
Integ(jj)=trapz(t,fun1);
end
K_f1=(1/sqrt(pi))*exp(j*(k*rho*alpha+pi/4)).*Integ*exp(-j*k*rho);
re=-1/n*cos(beta/2).*cot((pi+beta)/(2*n)).*K_f1;