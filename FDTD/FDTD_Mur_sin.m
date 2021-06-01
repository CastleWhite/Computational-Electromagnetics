%***********************************************************************
%     3-D FDTD code with Mur boundaries
%***********************************************************************
%
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a three-dimensional
%     Cartesian space lattice comprised of uniform cubic grid cells.
%     
%     The length, width, and height of the FDTD area are 
%     20.0 cm (x-direction), 20.0 cm (y-direction), and 
%     6.0 cm (z-direction), respectively.
%
%     The area is excited by an additive current source oriented
%     along the z-direction.  The source waveform is a sinusoid wave 
%     given by 
%          E(t)=ea*sin(2*pi*f).
%
%     The spatial step is 0.2cm.
%
%     This M-file displays the FDTD-computed Ez fields at every other
%     time step, and records those frames in a movie matrix, M, which 
%     is played at the end of the simulation using the "movie" command.
%
%***********************************************************************

clear

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=100;       %number of grid cells in x-direction
je=100;        %number of grid cells in y-direction
ke=30;        %number of grid cells in z-direction

ib=ie+1;     
jb=je+1;   
kb=ke+1;   

is=50;        %location of z-directed current source
js=50;        %location of z-directed current source

kobs=15;
dx=0.002;          %space increment of cubic lattice
dt=dx/(2.0*cc);    %time step

nmax=500;          %total number of time steps

%***********************************************************************
%     Sinusoid wave excitation
%***********************************************************************

f=2000e+6;
ea=-5e+5;

%***********************************************************************
%     Material parameters
%***********************************************************************

eps=1.0;
sig=0.0;        

%***********************************************************************
%     Updating coefficients
%***********************************************************************

ca=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
cb=(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));
da=1.0;
db=dt/muz/dx;

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(ie,jb,kb);
ey=zeros(ib,je,kb);
ez=zeros(ib,jb,ke);
hx=zeros(ib,je,ke);
hy=zeros(ie,jb,ke);
hz=zeros(ie,je,kb);

ex2=zeros(ie,jb,kb);   %the intermediate variable of Mur boundaries
ey2=zeros(ib,je,kb);   %the intermediate variable of Mur boundaries
ez2=zeros(ib,jb,ke);   %the intermediate variable of Mur boundaries


%***********************************************************************
%     Movie initialization
%***********************************************************************

tview(:,:)=ez(:,:,kobs);
sview(:,:)=ez(:,js,:);

subplot('position',[0.15 0.45 0.7 0.45]),pcolor(tview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['Ez(i,j,k=5), time step = 0']);
xlabel('i coordinate');
ylabel('j coordinate');

subplot('position',[0.15 0.10 0.7 0.25]),pcolor(sview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['Ez(i,j=13,k), time step = 0']);
xlabel('i coordinate');
ylabel('k coordinate');

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/2,gcf,rect);

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax
   
%***********************************************************************
%     Update electric fields
%***********************************************************************

ex(1:ie,2:je,2:ke)=ca*ex(1:ie,2:je,2:ke)+...
                   cb*(hz(1:ie,2:je,2:ke)-hz(1:ie,1:je-1,2:ke)+...
                       hy(1:ie,2:je,1:ke-1)-hy(1:ie,2:je,2:ke));

ey(2:ie,1:je,2:ke)=ca*ey(2:ie,1:je,2:ke)+...
                   cb*(hx(2:ie,1:je,2:ke)-hx(2:ie,1:je,1:ke-1)+...
                       hz(1:ie-1,1:je,2:ke)-hz(2:ie,1:je,2:ke));
                    
ez(2:ie,2:je,1:ke)=ca*ez(2:ie,2:je,1:ke)+...
                   cb*(hx(2:ie,1:je-1,1:ke)-hx(2:ie,2:je,1:ke)+...
                       hy(2:ie,2:je,1:ke)-hy(1:ie-1,2:je,1:ke));
                    
ez(is,js,1:ke)=ez(is,js,1:ke)+...
               ea*sin(2*pi*f*n*dt);
           
%**************************************************************************
%                          Mur boundaries
%**************************************************************************

ey(1,2:je-1,2:ke)=ey2(2,2:je-1,2:ke)+((cc*dt-dx)/(cc*dt+dx))*(ey(2,2:je-1,2:ke)-...
                  ey2(1,2:je-1,2:ke));   % surface x=0, Ey
         
ez(1,2:je,2:ke-1)=ez2(2,2:je,2:ke-1)+((cc*dt-dx)/(cc*dt+dx))*(ez(2,2:je,2:ke-1)-...
                  ez2(1,2:je,2:ke-1));   % surface x=0, Ez
              
ey(ib,2:je-1,2:ke)=ey2(ie,2:je-1,2:ke)+((cc*dt-dx)/(cc*dt+dx))*(ey(ie,2:je-1,2:ke)-...
                  ey2(1,2:je-1,2:ke));   % surface x=20.0cm, Ey
         
ez(ib,2:je,2:ke-1)=ez2(ie,2:je,2:ke-1)+((cc*dt-dx)/(cc*dt+dx))*(ez(ie,2:je,2:ke-1)-...
                  ez2(1,2:je,2:ke-1));   % surface x=20.0cm, Ez
              
ex(2:ie-1,1,2:ke)=ex2(2:ie-1,2,2:ke)+((cc*dt-dx)/(cc*dt+dx))*(ex(2:ie-1,2,2:ke)-...
                  ex2(2:ie-1,1,2:ke));   % surface y=0, Ex
         
ez(2:ie,1,2:ke-1)=ez2(2:ie,2,2:ke-1)+((cc*dt-dx)/(cc*dt+dx))*(ez(2:ie,2,2:ke-1)-...
                  ez2(2:ie,1,2:ke-1));   % surface y=0, Ez
              
ex(2:ie-1,jb,2:ke)=ex2(2:ie-1,je,2:ke)+((cc*dt-dx)/(cc*dt+dx))*(ex(2:ie-1,je,2:ke)-...
                  ex2(2:ie-1,jb,2:ke));   % surface y=10.0cm, Ex
         
ez(2:ie,jb,2:ke-1)=ez2(2:ie,je,2:ke-1)+((cc*dt-dx)/(cc*dt+dx))*(ez(2:ie,je,2:ke-1)-...
                  ez2(2:ie,jb,2:ke-1));   % surface y=10.0cm, Ez
              
ex(2:ie-1,2:je,1)=ex2(2:ie-1,2:je,2)+((cc*dt-dx)/(cc*dt+dx))*(ex(2:ie-1,2:je,2)-...
                  ex2(2:ie-1,2:je,1));    % surface z=0, Ex 
              
ey(2:ie,2:je-1,1)=ey2(2:ie,2:je-1,2)+((cc*dt-dx)/(cc*dt+dx))*(ey(2:ie,2:je-1,2)-...
                  ey2(2:ie,2:je-1,1));    % surface z=0, Ey
              
ex(2:ie-1,2:je,kb)=ex2(2:ie-1,2:je,ke)+((cc*dt-dx)/(cc*dt+dx))*(ex(2:ie-1,2:je,ke)-...
                  ex2(2:ie-1,2:je,kb));   % surface z=6.0cm, Ex

ey(2:ie,2:je-1,kb)=ey2(2:ie,2:je-1,ke)+((cc*dt-dx)/(cc*dt+dx))*(ey(2:ie,2:je-1,ke)-...
                  ey2(2:ie,2:je-1,kb));   % surface z=0.6cm, Ey
              
%***********************************************************************
%     Update magnetic fields
%***********************************************************************

hx(2:ie,1:je,1:ke)=hx(2:ie,1:je,1:ke)+...
                   db*(ey(2:ie,1:je,2:kb)-ey(2:ie,1:je,1:ke)+...
                       ez(2:ie,1:je,1:ke)-ez(2:ie,2:jb,1:ke));
                
hy(1:ie,2:je,1:ke)=hy(1:ie,2:je,1:ke)+...
                   db*(ex(1:ie,2:je,1:ke)-ex(1:ie,2:je,2:kb)+...
                       ez(2:ib,2:je,1:ke)-ez(1:ie,2:je,1:ke));
                
hz(1:ie,1:je,2:ke)=hz(1:ie,1:je,2:ke)+...
                   db*(ex(1:ie,2:jb,2:ke)-ex(1:ie,1:je,2:ke)+...
                       ey(1:ie,1:je,2:ke)-ey(2:ib,1:je,2:ke));
                    
%***********************************************************************
%     Visualize fields
%***********************************************************************

if mod(n,2)==0;

timestep=int2str(n);
tview(:,:)=ez(:,:,kobs);
sview(:,:)=ez(:,js,:);

subplot('position',[0.15 0.45 0.7 0.45]),pcolor(tview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['Ez(i,j,k=5), time step = ',timestep]);
xlabel('i coordinate');
ylabel('j coordinate');

subplot('position',[0.15 0.10 0.7 0.25]),pcolor(sview');
shading flat;
caxis([-1.0 1.0]);
colorbar;
axis image;
title(['Ez(i,j=13,k), time step = ',timestep]);
xlabel('i coordinate');
ylabel('k coordinate');

nn=n/2;
M(:,nn)=getframe(gcf,rect);

end;

%***********************************************************************
%                    Update intermediate variables
%***********************************************************************

ex2(1:ie,1:jb,1:kb)=ex(1:ie,1:jb,1:kb);
ey2(1:ib,1:je,1:kb)=ey(1:ib,1:je,1:kb);
ez2(1:ib,1:jb,1:ke)=ez(1:ib,1:jb,1:ke);

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

movie(gcf,M,0,10,rect);
