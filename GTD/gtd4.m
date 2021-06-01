clear all;
 
epsr    = 10; 
c       = 2.99792458e8;
R       = 0.01; 
n       = 1;
m       = sqrt(epsr);
 
order1  = n - 1/2;
order2  = n + 1/2;
 
N_re = 100;
N_im = 100;
f_re_min = 6.470e9;
f_re_max = 6.480e9;
f_im_min = 76.6e7;
f_im_max = 76.9e7;
f_re_step = (f_re_max - f_re_min)/(N_re - 1);
f_im_step = (f_im_max - f_im_min)/(N_im - 1);
 
for n_re = 1:N_re
    for n_im = 1:N_im
        f_re(n_re)  = f_re_min + (n_re-1)*f_re_step;
        f_im(n_im)  = f_im_min + (n_im-1)*f_im_step;
        f = f_re(n_re) + j*f_im(n_im);
        alpha       = 2*pi*R/(c/f);
        %   TE
        t1 = besselj(order1, m*alpha)/besselj(order2, m*alpha);
        t2 = besselh(order1, 2, alpha)/besselh(order2, 2, alpha)*m;
        errorTM = t1 - t2 + (n/alpha)*(m^2-1)/m;
        DET(n_re, n_im) = errorTM;        
    end
end
    
%   Find resonant frequency and Q factor
[row,I] = min(DET);
[MIN,J] = min(row);
REAL_F  = f_re(I(J))
IMAG_F  = f_im(J)
Q       = REAL_F/(2*IMAG_F)
mesh(f_im, f_re, 1./abs(DET));
