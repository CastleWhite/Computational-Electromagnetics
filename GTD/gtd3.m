clear;clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
c=3*10^8; % 波速 
f=3*10^9; % 频率 
lamda=c/f; % 波长 
k=2*pi/lamda; % 波数 
 
epsz=1/(4*pi*9*10^9); % 真空介电常数 
mu=4*pi*10.^(-7); % 真空磁导率 
Z=120*pi; % 真空波阻抗 
epsilon=1; % 相对介电常数 
sigma=0; % 电导率 
 
N=100; % 网格数量 
L=800; % 迭代次数 
ddx=lamda/20; % 网格尺寸 
dt=ddx/(2*c); % 时间间隔 
 
ia=N/4; % 总场区域x左 
ib=3*N/4; % 总场区域x右 
ja=ia; % 总场区域x下 
jb=ib; % 总场区域x上 
 
npml=N/8; % PML点数 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
r=2*lamda;  
M=r/(2*ddx); 
 
for i=1:N  
for j=1:N 
ga(i,j)=1/(epsilon+sigma*dt/epsz); % 求和参量 
gb(i,j)=sigma*dt/epsz; % 求和参量 
end 
end 
 
spread=8; % 脉冲宽度 
t0=25; % 脉冲高度 
is=N/2; % 源的X位置 
js=N/2; % 源的Y位置 
 
ez_inc=zeros(1,N); 
hx_inc=zeros(1,N); 
 
ez_low1=0; 
ez_low2=0; 
ez_low3=0; 
ez_low4=0; 
 
dz=zeros(N,N); % z方向电荷密度 
ez=zeros(N,N); % z方向电场 
iz=zeros(N,N); % z方向电场求和参量 
hx=zeros(N,N); % x方向磁场 
hy=zeros(N,N); % y方向磁场 
ihx=zeros(N,N); % x方向磁场参量 
ihy=zeros(N,N); % y方向磁场参量 
 
%%%%%%%%%%%%%%%%%%%%%%%%PML%%%%%%%%%%%%%%%%%%% 
for i=1:N 
gi2(i)=1; 
gi3(i)=1; 
fi1(i)=0; 
fi2(i)=1; 
fi3(i)=1; 
end 
for j=1:N 
gj2(j)=1; 
gj3(j)=1; 
fj1(j)=0; 
fj2(j)=1; 
fj3(j)=1; 
end 
 
for i=1:npml+1 
xnum=npml-i+1; 
xxn=xnum/npml;  
xn=0.33*(xxn^3); 
gi2(i)=1/(1+xn); 
gi2(N-i+1)=1/(1+xn); 
gi3(i)=(1-xn)/(1+xn); 
gi3(N-i+1)=(1-xn)/(1+xn); 
 
xxn=(xnum-0.5)/npml; 
xn=0.25*(xxn^3); 
fi1(i)=xn; 
fi1(N-i)=xn; 
fi2(i)=1/(1+xn); 
fi2(N-i)=1/(1+xn); 
fi3(i)=(1-xn)/(1+xn); 
fi3(N-i)=(1-xn)/(1+xn); 
end 
 
for j=1:npml+1 
xnum=npml-j+1; 
xxn=xnum/npml;  
xn=0.33*(xxn^3); 
gj2(j)=1/(1+xn); 
gj2(N-j+1)=1/(1+xn); 
gj3(j)=(1-xn)/(1+xn); 
gj3(N-j+1)=(1-xn)/(1+xn); 
 
xxn=(xnum-0.5)/npml; 
xn=0.25*(xxn^3); 
fj1(j)=xn; 
fj1(N-j)=xn; 
fj2(j)=1/(1+xn); 
fj2(N-j)=1/(1+xn); 
fj3(j)=(1-xn)/(1+xn); 
fj3(N-j)=(1-xn)/(1+xn); 
end 
 
for T=1:L 
%%%%%%%%%%%%%%TM波Y方向传播%%%%%%%%%%%%%%%%%%%%% 
for j=2:N-1 
ez_inc(j)=ez_inc(j)+0.5*(hx_inc(j-1)-hx_inc(j)); 
end 
 
ez_inc(1)=ez_low2; 
ez_low2=ez_low1; 
ez_low1=ez_inc(2); 
ez_inc(N)=ez_low3; 
ez_low3=ez_low4; 
ez_low4=ez_inc(N-1); 
%% %%%%%%%%%%%%%% 电荷密度dz%%%%%%%%%%%%%%%%%%%%%%%%% 
for i=2:N-1  
for j=2:N 
dz(i,j)=gi3(i)*gj3(j)*dz(i,j)+gi2(i)*gj2(j)*0.5*( hy(i,j)-hy(i-1,j)-hx(i,j)+hx(i,j-1) ); 
end 
end  
%%%%%%%%%%%%%%%%%%脉冲的加入%%%%%%%%%%%%%%%%%%% 
source=exp((-0.5)*( (t0-T)/spread ).^2);  
ez_inc(5)=source;  
 
for i=ia:ib 
dz(i,ja)=dz(i,ja)+0.5*hx_inc(ja-1); 
dz(i,jb)=dz(i,jb)-0.5*hx_inc(jb); 
end 
%%%%%%%%%%%%%% 电场ez%%%%%%%%%%%%%%%%%%%%%%%% 
for i=1:N  
for j=1:N 
ez(i,j)=ga(i,j)*( dz(i,j)-iz(i,j) ); 
iz(i,j)=iz(i,j)+gb(i,j)*ez(i,j) ; 
end 
end  
%%%%%%%%%%%%%%%%%%%%边界电场%%%%%%%%%%%%%%%%%%%%%%% 
for j=1:N 
ez(1,j)=0; 
ez(N,j)=0; 
end 
 
for i=1:N 
ez(i,1)=0; 
ez(i,N)=0; 
end  
 
%%%%%%%%%%%%%%%%%%%%%%边界条件%%%%%%%%%%%%%%%%%%%%%%%% 
for i=N/2-M:N/2+M-1 
for j=N/2-M:N/2+M-1 
ez(i,j)=0; 
end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%磁场%%%%%%%%%%%%%%%%%%%%%%%% 
for j=1:N-1 
hx_inc(j)=hx_inc(j)+0.5*(ez_inc(j)-ez_inc(j+1)); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%磁场Hx%%%%%%%%%%%%%%%%%%%%%% 
for i=1:N  
for j=1:N-1 
curl_e=ez(i,j)-ez(i,j+1); 
ihx(i,j)=ihx(i,j)+fi1(i)*curl_e; 
hx(i,j)=fj3(j)*hx(i,j)+fj2(j)*0.5*(curl_e+ihx(i,j)); 
end 
end;  
 
for i=ia:ib 
hx(i,ja-1)=hx(i,ja-1)+0.5*ez_inc(ja); 
hx(i,jb)=hx(i,jb)-0.5*ez_inc(jb); 
end 
%%%%%%%%磁场Hy%%%%%%%%%%%%%%%%%%% 
for i=1:N-1  
for j=1:N 
curl_e=ez(i+1,j)-ez(i,j); 
ihy(i,j)=ihy(i,j)+fj1(j)*curl_e; 
hy(i,j)=fi3(i)*hy(i,j)+fi2(i)*0.5*(curl_e+ihy(i,j)); 
end 
end  
 
for j=ja:jb; 
hy(ia-1,j)=hy(ia-1,j)-0.5*ez_inc(j); 
hy(ib,j)=hy(ib,j)+0.5*ez_inc(j); 
end 
 
mesh(ez) 
drawnow; 
end
