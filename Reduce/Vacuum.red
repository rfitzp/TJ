% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calulation of TJ vacuum solution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out "Vacuum.out"$

% %%%%%%%%%%%%%
linelength 120$
% %%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to average var**2 terms in general polynomial ex by
% performing integral over general angular variable ivar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure av2(ex,var,ivar)$
   begin scalar c0,c1,c2$
   c0 := coeffn(ex,var,0)$
   c1 := coeffn(ex,var,1)$
   c2 := coeffn(ex,var,2)$
   c2 := int(c2,ivar)$
   return c0 + var*c1 + var**2*(sub(ivar=2*pi,c2) - sub(ivar=0,c2))/(2*pi)$
end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express inverse of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure inv(ex,var)$
   begin scalar c0,c1,c2,i0,i1,i2$
   c0 := coeffn(ex,var,0)$
   c1 := coeffn(ex,var,1)/c0$
   c2 := coeffn(ex,var,2)/c0$
   i0 := 1/c0$
   i1 := -c1/c0$
   i2 := (c1**2 - c2)/c0$
   return i0 + var*i1 + var**2*i2$
end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to return amplitude of nth-order Fourier cosine
% component in angle-like variable ivar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure fcos(ex,ivar,n)$
begin scalar x,y,norm$
   x := ex*cos(n*ivar)$
   y := int(x,ivar)$
   if n=0 then norm := 1 else norm := 2$
   return norm*(sub(ivar=2*pi,y) - sub(ivar=0,y))/(2*pi)
end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to return amplitude of nth-order Fourier sine
% component in angle-like variable ivar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure fsin(ex,ivar,n)$
   begin scalar x,y,norm$
   x := ex*sin(n*ivar)$
   y := int(x,ivar)$
   if n=0 then norm := 1 else norm := 2$
   return norm*(sub(ivar=2*pi,y) - sub(ivar=0,y))/(2*pi)
   end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express nth power of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure pow1(ex,var,n)$
   begin scalar c0,c1,c2,i0,i1,i2$
      c0 := coeffn(ex,var,0)$
      c1 := coeffn(ex,var,1)/c0$
      c2 := coeffn(ex,var,2)/c0$
      i1 := n*c1$
      i2 := n*(n-1)*c1*c1/2 + n*c2$
      return c0**n*(1 + var*i1 + var**2*i2)
   end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express n/m th power of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure pow2(ex,var,n,m)$
   begin scalar c0,c1,c2,i0,i1,i2$
      c0 := coeffn(ex,var,0)$
      c1 := coeffn(ex,var,1)/c0$
      c2 := coeffn(ex,var,2)/c0$
      i1 := n*c1/m$
      i2 := n*(n-m)*c1*c1/2/m/m + n*c2/m$
      return c0**(n/m)*(1 + var*i1 + var**2*i2)
   end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express square root of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure sqr(ex,var)$
   begin scalar c0,c1,c2,i0,i1,i2$
      c0 := coeffn(ex,var,0)$
      c1 := coeffn(ex,var,1)/c0$
      c2 := coeffn(ex,var,2)/c0$
      i1 := c1/2$
      i2 := - c1*c1/8 + c2/2$
      return c0**(1/2)*(1 + var*i1 + var**2*i2)
end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful trigonometrical identities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for all x, y let cos(x) * cos(y) = ( cos(x+y) + cos(x-y) )/2,
    	     	 cos(x) * sin(y) = ( sin(x+y) - sin(x-y) )/2,
		 sin(x) * sin(y) = ( cos(x-y) - cos(x+y) )/2$
for all x    let sin(x)**2 = ( 1 - cos(2*x) )/2,
                 cos(x)**2 = ( 1 + cos(2*x) )/2$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that inverse aspect ratio eps is ordered first, both
% internally and externally
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
korder eps,i$
order  eps,i$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that all quantities expressed as polynomials in eps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
factor eps, cos, sin$
on  rat$
on  revpri$
off ratpri$

% %%%%%%%%%%%%%%%%%%%%%%%%
% Truncate at order eps**2
% %%%%%%%%%%%%%%%%%%%%%%%%
let eps**3 = 0$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define equilibirum quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depend H1,r$
depend H2,r$
depend V2,r$
depend H3,r$
depend V3,r$
depend H4,r$
depend V4,r$
depend lg,r$

H1pp := df(H1,r)/r - 1$
H2pp := df(H2,r)/r + (2*2-1)*H2/r**2$
H3pp := df(H3,r)/r + (3*3-1)*H3/r**2$
H4pp := df(H4,r)/r + (4*4-1)*H4/r**2$
V2pp := df(V2,r)/r + (2*2-1)*V2/r**2$
V3pp := df(V3,r)/r + (3*3-1)*V3/r**2$
V4pp := df(V4,r)/r + (4*4-1)*V4/r**2$

gr21 := + 2*df(H1,r)*cos(t)
        + 2*df(H2,r)*cos(2*t)
	+ 2*df(H3,r)*cos(3*t)
	+ 2*df(H4,r)*cos(4*t)
        + 2*df(V2,r)*sin(2*t)
	+ 2*df(V3,r)*sin(3*t)
	+ 2*df(V4,r)*sin(4*t)$
gr22 :=    3*r**2/4 - H1
	+ (df(H1,r)*df(H1,r) + (1*1-1)*H1*H1/r/r)/2
	+ (df(H2,r)*df(H2,r) + (2*2-1)*H2*H2/r/r)/2
	+ (df(H3,r)*df(H3,r) + (3*3-1)*H3*H3/r/r)/2
	+ (df(H4,r)*df(H4,r) + (4*4-1)*H4*H4/r/r)/2
	+ (df(V2,r)*df(V2,r) + (2*2-1)*V2*V2/r/r)/2
	+ (df(V3,r)*df(V3,r) + (3*3-1)*V3*V3/r/r)/2
	+ (df(V4,r)*df(V4,r) + (4*4-1)*V4*V4/r/r)/2$

grr := 1 + eps*gr21 + eps**2*gr22$

grt21 :=   sin(t)  
         - (H1pp + df(H1,r)/r + (1*1-1)*H1/r/r)*sin(t)  /1             
	 - (H2pp + df(H2,r)/r + (2*2-1)*H2/r/r)*sin(2*t)/2
         - (H3pp + df(H3,r)/r + (3*3-1)*H3/r/r)*sin(3*t)/3
         - (H4pp + df(H4,r)/r + (4*4-1)*H4/r/r)*sin(4*t)/4
	 + (V2pp + df(V2,r)/r + (2*2-1)*V2/r/r)*cos(2*t)/2
	 + (V3pp + df(V3,r)/r + (3*3-1)*V3/r/r)*cos(3*t)/3
	 + (V4pp + df(V4,r)/r + (4*4-1)*V4/r/r)*cos(4*t)/4$
grt22 :=   (2-1)*df(V2,r)*H2/2/r/r - (2-1)*df(H2,r)*V2/2/r/r
         + (3-1)*df(V3,r)*H3/2/r/r - (3-1)*df(H3,r)*V3/2/r/r
	 + (4-1)*df(V4,r)*H4/2/r/r - (4-1)*df(H4,r)*V4/2/r/r
	 + (2-1)*V2pp    *H2/2/2/r - (2-1)         *H2pp*V2/2/2/r 
	 + (3-1)*V3pp    *H3/2/3/r - (3-1)         *H3pp*V3/2/3/r
	 + (4-1)*V4pp    *H4/2/4/r - (4-1)         *H4pp*V4/2/4/r
	 +       V2pp*df(H2,r)/2/2 - H2pp*df(V2,r)/2/2 
	 +       V3pp*df(H3,r)/2/3 - H3pp*df(V3,r)/2/3
	 +       V4pp*df(H4,r)/2/4 - H4pp*df(V4,r)/2/4$

grt := eps*grt21 + eps*eps*grt22$

R21 := - 2*r*cos(t)$
R22 := - (r*r/2 - 2*H1 - r*df(H1,r))$

R2 := 1 + eps*R21 + eps**2*R22$

hfun := av2(grr*R2,eps,t)$
ifun := i*r*av2(grt*R2,eps,t)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate transformation of angles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fw := r*sin(t)
   -       (df(H1,r) - (1-1)*H1/r) * sin(t)
   - (1/2)*(df(H2,r) - (2-1)*H2/r) * sin(2*t)
   - (1/3)*(df(H3,r) - (3-1)*H3/r) * sin(3*t)
   - (1/4)*(df(H4,r) - (4-1)*H4/r) * sin(4*t)
   + (1/2)*(df(V2,r) - (2-1)*V2/r) * cos(2*t)
   + (1/3)*(df(V3,r) - (3-1)*V3/r) * cos(3*t)
   + (1/4)*(df(V4,r) - (4-1)*V4/r) * cos(4*t)$

Feta := - (r/2)*sin(eta)
   + cos(pi)  *(H1/r)*sin(eta)
   + cos(2*pi)*(H2/r)*sin(2*eta)
   + cos(3*pi)*(H3/r)*sin(3*eta)
   + cos(4*pi)*(H4/r)*sin(4*eta)
   + cos(2*pi)*(V2/r)*cos(2*eta)
   + cos(3*pi)*(V3/r)*cos(3*eta)
   + cos(4*pi)*(V4/r)*cos(4*eta)$

Ftheta := - Fw - sub(eta=pi-t,Feta)$

Ftheta_target := - (r/2)*sin(t)
   +       (df(H1,r) + H1/r) * sin(t)
   + (1/2)*(df(H2,r) + H2/r) * sin(2*t)
   + (1/3)*(df(H3,r) + H3/r) * sin(3*t)
   + (1/4)*(df(H4,r) + H4/r) * sin(4*t)
   - (1/2)*(df(V2,r) + V2/r) * cos(2*t)
   - (1/3)*(df(V3,r) + V3/r) * cos(3*t)
   - (1/4)*(df(V4,r) + V4/r) * cos(4*t)$

write ("Ftheta residual:");
Ftheta - Ftheta_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate toroidal angle eta as function of straight field-line angle t
%
%  eta = pi - t - eps*Ftheta(t)
%
% Express all trigonometric functions in terms of t
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt := pi-t$
let           cos(eta)   = cos(tt)   +     sin(tt)*eps*Ftheta$
for all n let cos(n*eta) = cos(n*tt) + n*sin(n*tt)*eps*Ftheta$
let           sin(eta)   = sin(tt)   -     cos(tt)*eps*Ftheta$
for all n let sin(n*eta) = sin(n*tt) - n*cos(n*tt)*eps*Ftheta$

% %%%%%%%%%%%%
% Calculate xi
% %%%%%%%%%%%%
xi := 1 + eps*(
   (r/2)*cos(eta)
   + cos(pi)  *(H1/r)*cos(eta)
   + cos(2*pi)*(H2/r)*cos(2*eta) - cos(2*pi)*(V2/r)*sin(2*eta)
   + cos(3*pi)*(H3/r)*cos(3*eta) - cos(3*pi)*(V3/r)*sin(3*eta)
   + cos(4*pi)*(H4/r)*cos(4*eta) - cos(4*pi)*(V4/r)*sin(4*eta))
   + eps*eps*(5*r**2/16 - H1/4
   + (3/4)*((H1/r)**2 + (H2/r)**2 + (H3/r)**2 + (H4/r)**2)
   + (3/4)*(            (V2/r)**2 + (V3/r)**2 + (V4/r)**2))$
xi  := av2(xi,eps,t)$
xi0 := coeffn(xi,eps,0)$
xi1 := coeffn(xi,eps,1)$
xi2 := coeffn(xi,eps,2)$

xi0_target := 1$
xi1_target := - (r/2)*cos(t)
   + (H1/r)*cos(t)
   + (H2/r)*cos(2*t) + (V2/r)*sin(2*t)
   + (H3/r)*cos(3*t) + (V3/r)*sin(3*t)
   + (H4/r)*cos(4*t) + (V4/r)*sin(4*t)$
xi2_target := (3/16)*r**2 + H1/4 + r*df(H1,r)/4
   + (1/4)*((H1/r)**2 + (H2/r)**2 + (H3/r)**2 + (H4/r)**2)
   + (1/4)*(            (V2/r)**2 + (V3/r)**2 + (V4/r)**2)
   - (1/2)*(df(H1,r)*H1/r + df(H2,r)*H2/r + df(H3,r)*H3/r + df(H4,r)*H4/r)
   - (1/2)*(df(V1,r)*V1/r + df(V2,r)*V2/r + df(V3,r)*V3/r + df(V4,r)*V4/r)$

write ("xi0 residual:");
xi0 - xi0_target;
write ("xi1 residual:");
xi1 - xi1_target;
write ("xi2 residual:");
xi2 - xi2_target;

% %%%%%%%%%%%%%%
% Calculate lnxi
% %%%%%%%%%%%%%%
lnxi0 := 0$
lnxi1 := xi1$
lnxi2 := xi2 - av2(xi1*xi1,eps,t)/2$
lnxi  := eps*lnxi1 + eps**2*lnxi2$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate vacuum coupling coefficients for Z
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let s*s = 1$

xir  := xi/r$
rxi  := r*inv(xi,eps)$

f1    := 1 - eps*rxi*cos(eta)$
fac1  := sqr(f1,eps)$
fac2  := pow2(xi,eps,m,1)$
fac2a := inv(fac2,eps)$
fac3  := 1 - eps*eps*((n-m+1/2)*(n-m+3/2)/4/(m-1) + n/2)*r**2$
fac3a := 1 + eps*eps*((m-n+1/2)*(m-n+3/2)/4/(m+1) + n/2)*r**2$
fac4  := (1 - eps**2*m**2*Ftheta**2/2)*(cos(kk*t) - i*s*sin(kk*t))
   + eps*m*Ftheta*(sin(kk*t) + i*s*cos(kk*t))$

fac0  := lg + lnxi + eps**2*(1/4)*(n**2-5/4 + (n**2+3/4)*lg)*r**2$
fac01 := cos(kk*t) - i*sin(kk*t)$

fac10 := 1 - eps*eps*((n**2-1/4)*lg/2 + (n**2+3/4)/4)*r**2$
fac11 := (1 - eps**2*Ftheta**2/2)*(cos(kk*t) - i*s*sin(kk*t))
   + eps*Ftheta*(sin(kk*t) + i*s*cos(kk*t))$

Pfac := fac1*fac2 *fac3 *fac4$
Qfac := fac1*fac2a*fac3a*fac4$

Pfac0 := fac1*fac0*fac01$
Pfac1 := fac1*xi*fac10*fac11$

hfac  := hfun * (cos(kk*t) + i*s*sin(kk*t))$
ifac  := ifun * (cos(kk*t) + i*s*sin(kk*t))$

off   revpri$
order s,m,n,r$
order H1,df(H1,r),df(H2,r,2),H2,df(V2,r,2),V2,df(H3,r,2),H3,df(V3,r,2),V3,df(H4,r,2),H4,df(V4,r,2),V4$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=0 coupling coefficients 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = 0$
Pmm := Pfac$
Qmm := Qfac$
P00 := Pfac0$
P01 := Pfac1$
hh0 := hfac$
ii0 := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmm0 := coeffn(Pmm,eps,0)$
Qmm0 := coeffn(Qmm,eps,0)$
Pmm2 := coeffn(Pmm,eps,2)$
Qmm2 := coeffn(Qmm,eps,2)$

Pmm0_target := 1$
Qmm0_target := 1$
Pmm2_target :=
          i*s*m**2*
                (  (df(H2,r)*V2 - df(V2,r)*H2)/2/2/r
      	         + (df(H3,r)*V3 - df(V3,r)*H3)/2/3/r
	         + (df(H4,r)*V4 - df(V4,r)*H4)/2/4/r)
        - m**2* (  (df(H1,r)**2              )/(4*1**2)
	         + (df(H2,r)**2 + df(V2,r)**2)/(4*2**2)
	         + (df(H3,r)**2 + df(V3,r)**2)/(4*3**2)
	         + (df(H4,r)**2 + df(V4,r)**2)/(4*4**2))
        - m**2* (  (df(H1,r)*H1/r                )/(2*1**2)
	         + (df(H2,r)*H2/r + df(V2,r)*V2/r)/(2*2**2)
		 + (df(H3,r)*H3/r + df(V3,r)*V3/r)/(2*3**2)
                 + (df(H4,r)*H4/r + df(V4,r)*V4/r)/(2*4**2))
        + m**2* (  (2*2-1)* ((H2/r)**2 + (V2/r)**2)/(4*2**2)
	         + (3*3-1)* ((H3/r)**2 + (V3/r)**2)/(4*3**2)
	         + (4*4-1)* ((H4/r)**2 + (V4/r)**2)/(4*4**2))
        - m*    (  (df(H1,r)*H1/r                )/2
		 + (df(H2,r)*H2/r + df(V2,r)*V2/r)/2
		 + (df(H3,r)*H3/r + df(V3,r)*V3/r)/2
		 + (df(H4,r)*H4/r + df(V4,r)*V4/r)/2)
	+       (m**2/4 + (m-1)/4)*r*df(H1,r)
        +       (3*m/4 - 1/2)*H1
        +       (- 4*m**2 + 11*m - 4*n**2 - 6)*r**2/16/(m-1)$

G2 :=  - r*df(H1,r)/4
       + (df(H1,r)**2 + 2*df(H1,r)*H1/r - (1*1-1)*(H1/r)**2)/(4*1**2)
       + (df(H2,r)**2 + 2*df(H2,r)*H2/r - (2*2-1)*(H2/r)**2)/(4*2**2)
       + (df(V2,r)**2 + 2*df(V2,r)*V2/r - (2*2-1)*(V2/r)**2)/(4*2**2)
       + (df(H3,r)**2 + 2*df(H3,r)*H3/r - (3*3-1)*(H3/r)**2)/(4*3**2)
       + (df(V3,r)**2 + 2*df(V3,r)*V3/r - (3*3-1)*(V3/r)**2)/(4*3**2)
       + (df(H4,r)**2 + 2*df(H4,r)*H4/r - (4*4-1)*(H4/r)**2)/(4*4**2)
       + (df(V4,r)**2 + 2*df(V4,r)*V4/r - (4*4-1)*(V4/r)**2)/(4*4**2)
       - i*s*(  (df(H2,r)*V2/r - df(V2,r)*H2/r)/2
      	      + (df(H3,r)*V3/r - df(V3,r)*H3/r)/3
	      + (df(H4,r)*V4/r - df(V4,r)*H4/r)/4)/2$
G1 :=  - 3*H1/4 - r*df(H1,r)/4
       + (df(H1,r)*H1/r                )/2
       + (df(H2,r)*H2/r + df(V2,r)*V2/r)/2
       + (df(H3,r)*H3/r + df(V3,r)*V3/r)/2
       + (df(H4,r)*H4/r + df(V4,r)*V4/r)/2$
G0 := - ((n**2 + (m-2)*(m-3/4))/4/(m-1))*r**2 - H1/2 - r*df(H1,r)/4$

Pmm2_target := G0 - m*G1 - m*m*G2$

Qmm2_target := 
        - i*s*m**2*
                (  (df(H2,r)*V2 - df(V2,r)*H2)/2/2/r
      	         + (df(H3,r)*V3 - df(V3,r)*H3)/2/3/r
	         + (df(H4,r)*V4 - df(V4,r)*H4)/2/4/r)
        - m**2* (  (df(H1,r)**2              )/(4*1**2)
	         + (df(H2,r)**2 + df(V2,r)**2)/(4*2**2)
	         + (df(H3,r)**2 + df(V3,r)**2)/(4*3**2)
	         + (df(H4,r)**2 + df(V4,r)**2)/(4*4**2))
        - m**2* ( (df(H1,r)*H1/r                )/(2*1**2)
	         + (df(H2,r)*H2/r + df(V2,r)*V2/r)/(2*2**2)
		 + (df(H3,r)*H3/r + df(V3,r)*V3/r)/(2*3**2)
                 + (df(H4,r)*H4/r + df(V4,r)*V4/r)/(2*4**2))
        + m**2* (  (2*2-1)* ((H2/r)**2 + (V2/r)**2)/(4*2**2)
	         + (3*3-1)* ((H3/r)**2 + (V3/r)**2)/(4*3**2)
	         + (4*4-1)* ((H4/r)**2 + (V4/r)**2)/(4*4**2))
        + m*    (  (df(H1,r)*H1/r                )/2
		 + (df(H2,r)*H2/r + df(V2,r)*V2/r)/2
		 + (df(H3,r)*H3/r + df(V3,r)*V3/r)/2
		 + (df(H4,r)*H4/r + df(V4,r)*V4/r)/2)
	+       (m**2/4 + (-m-1)/4)*r*df(H1,r)
	+       (- 3*m/4 - 1/2)*H1
        +       (- 4*m**2 - 11*m - 4*n**2 - 6)*r**2/16/(-m-1)$

G2a := - r*df(H1,r)/4
       + (df(H1,r)**2 + 2*df(H1,r)*H1/r - (1*1-1)*(H1/r)**2)/(4*1**2)
       + (df(H2,r)**2 + 2*df(H2,r)*H2/r - (2*2-1)*(H2/r)**2)/(4*2**2)
       + (df(V2,r)**2 + 2*df(V2,r)*V2/r - (2*2-1)*(V2/r)**2)/(4*2**2)
       + (df(H3,r)**2 + 2*df(H3,r)*H3/r - (3*3-1)*(H3/r)**2)/(4*3**2)
       + (df(V3,r)**2 + 2*df(V3,r)*V3/r - (3*3-1)*(V3/r)**2)/(4*3**2)
       + (df(H4,r)**2 + 2*df(H4,r)*H4/r - (4*4-1)*(H4/r)**2)/(4*4**2)
       + (df(V4,r)**2 + 2*df(V4,r)*V4/r - (4*4-1)*(V4/r)**2)/(4*4**2)
       + i*s*(  (df(H2,r)*V2/r - df(V2,r)*H2/r)/2
      	      + (df(H3,r)*V3/r - df(V3,r)*H3/r)/3
	      + (df(H4,r)*V4/r - df(V4,r)*H4/r)/4)/2$
G0a := ((n**2 + (m+2)*(m+3/4))/4/(m+1))*r**2 - H1/2 - r*df(H1,r)/4$

Qmm2_target := G0a + m*G1 - m*m*G2a$

write ("Pmm0 residual:");
Pmm0 - Pmm0_target;
write ("Pmm2 residual:");
Pmm2 - Pmm2_target;
write ("Qmm0 residual:");
Qmm0 - Qmm0_target;
write ("Qmm2 residual:");
Qmm2 - Qmm2_target;

P000 := coeffn(P00,eps,0)$
P002 := coeffn(P00,eps,2)$

P000_target := lg$
P002_target := (n**2/4 - 5/16 + (n**2/4 + 3/8)*lg)*r**2 - (H1/2 + r*df(H1,r)/4)*lg - G1$

write ("P000 residual:");
P000 - P000_target;
write ("P002 residual:");
P002 - P002_target;

P010 := coeffn(P01,eps,0)$
P012 := coeffn(P01,eps,2)$

G0b         := - ((n**2/2 - 1/8)*lg + n**2/4)*r**2 - H1/2 - r*df(H1,r)/4$
P010_target := 1$
P012_target := G0b - G1 - G2$

write ("P010 residual:");
P010 - P010_target;
write ("P012 residual:");
P012 - P012_target;

hh00 := coeffn(hh0,eps,0)$
hh02 := coeffn(hh0,eps,2)$
ii00 := coeffn(ii0,eps,0)$
ii02 := coeffn(ii0,eps,2)$

hh00_target := 1$
hh02_target := r**2/4 + H1 - r*df(H1,r)
   + df(H1,r)**2/2
   + df(H2,r)**2/2 + df(V2,r)**2/2
   + df(H3,r)**2/2 + df(V3,r)**2/2
   + df(H4,r)**2/2 + df(V4,r)**2/2
   + (2*2-1)*(H2/r)**2/2 + (2*2-1)*(V2/r)**2/2
   + (3*3-1)*(H3/r)**2/2 + (3*3-1)*(V3/r)**2/2
   + (4*4-1)*(H4/r)**2/2 + (4*4-1)*(V4/r)**2/2$
ii00_target := 0$
ii02_target := 0$

write ("hh00 residual:");
hh00 - hh00_target;
write ("hh02 residual:");
hh02 - hh02_target;
write ("ii00 residual:");
ii00 - ii00_target;
write ("ii02 residual:");
ii02 - ii02_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = 1$
Pmmp1 := Pfac$
Qmmp1 := Qfac$
P00p1 := Pfac0$
P01p1 := Pfac1$
hhp1  := hfac$
iip1  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmp11 := coeffn(Pmmp1,eps,1)$
Qmmp11 := coeffn(Qmmp1,eps,1)$

Pmmp11_target := m * ((1/4/m - 1/2)*r + df(H1,r)/2 + H1/r)$
Qmmp11_target := m * ( r/4/m + df(H1,r)/2)$

write ("Pmmp11 residual:");
Pmmp11 - Pmmp11_target;
write ("Qmmp11 residual:");
Qmmp11 - Qmmp11_target;

P00p11 := coeffn(P00p1,eps,1)$

P00p11_target := lg*r/4 - r/4 + H1/2/r$

write ("P00p11 residual:");
P00p11 - P00p11_target;

P01p11 := coeffn(P01p1,eps,1)$

P01p11_target := ((1/4 - 1/2)*r + df(H1,r)/2 + H1/r)$

write ("P01p11 residual:");
P01p11 - P01p11_target;

hhp11 := coeffn(hhp1,eps,1)$
iip11 := coeffn(iip1,eps,1)$

hhp11_target := - r + df(H1,r)$
iip11_target := s*(-r + df(H1,r))$

write ("hhp11 residual:");
hhp11 - hhp11_target;
write ("iip11 residual:");
iip11 - iip11_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = 2$
Pmmp2 := Pfac$
Qmmp2 := Qfac$
P00p2 := Pfac0$
P01p2 := Pfac1$
hhp2  := hfac$
iip2  := ifac$

for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmp21 := coeffn(Pmmp2,eps,1)$
Qmmp21 := coeffn(Qmmp2,eps,1)$

Pmmp21_target := m * ((df(H2,r) - i*s*df(V2,r)) + (2+1)*(H2 - i*s*V2)/r)/(2*2)$
Qmmp21_target := m * ((df(H2,r) - i*s*df(V2,r)) - (2-1)*(H2 - i*s*V2)/r)/(2*2)$

write ("Pmmp21 residual:");
Pmmp21 - Pmmp21_target;
write ("Qmmp21 residual:");
Qmmp21 - Qmmp21_target;

P00p21 := coeffn(P00p2,eps,1)$

P00p21_target := (H2 - i*V2)/2/r$

write ("P00p21 residual:");
P00p21 - P00p21_target;

P01p21 := coeffn(P01p2,eps,1)$

P01p21_target := ((df(H2,r) - i*s*df(V2,r)) + (2+1)*(H2 - i*s*V2)/r)/(2*2)$

write ("P01p21 residual:");
P01p21 - P01p21_target;

hhp21 := coeffn(hhp2,eps,1)$
iip21 := coeffn(iip2,eps,1)$

hhp21_target :=    df(H2,r) + i*s*df(V2,r)$
iip21_target := s*(df(H2,r) + i*s*df(V2,r))/2 + s*(2*2-1)*(H2 + i*s*V2)/2/r$

write ("hhp21 residual:");
hhp21 - hhp21_target;
write ("iip21 residual:");
iip21 - iip21_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = 3$
Pmmp3 := Pfac$
Qmmp3 := Qfac$
P00p3 := Pfac0$
P01p3 := Pfac1$
hhp3  := hfac$
iip3  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmp31 := coeffn(Pmmp3,eps,1)$
Qmmp31 := coeffn(Qmmp3,eps,1)$

Pmmp31_target := m * ((df(H3,r) - i*s*df(V3,r)) + (3+1)*(H3 - i*s*V3)/r)/(2*3)$
Qmmp31_target := m * ((df(H3,r) - i*s*df(V3,r)) - (3-1)*(H3 - i*s*V3)/r)/(2*3)$

write ("Pmmp31 residual:");
Pmmp31 - Pmmp31_target;
write ("Qmmp31 residual:");
Qmmp31 - Qmmp31_target;

P00p31 := coeffn(P00p3,eps,1)$

P00p31_target := (H3 - i*V3)/2/r$

write ("P00p31 residual:");
P00p31 - P00p31_target;

P01p31 := coeffn(P01p3,eps,1)$

P01p31_target := ((df(H3,r) - i*s*df(V3,r)) + (3+1)*(H3 - i*s*V3)/r)/(2*3)$

write ("P01p31 residual:");
P01p31 - P01p31_target;

hhp31 := coeffn(hhp3,eps,1)$
iip31 := coeffn(iip3,eps,1)$

hhp31_target :=    df(H3,r) + i*s*df(V3,r)$
iip31_target := s*(df(H3,r) + i*s*df(V3,r))/3 + s*(3*3-1)*(H3 + i*s*V3)/3/r$

write ("hhp31 residual:");
hhp31 - hhp31_target;
write ("iip31 residual:");
iip31 - iip31_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = 4$
Pmmp4 := Pfac$
Qmmp4 := Qfac$
P00p4 := Pfac0$
P01p4 := Pfac1$
hhp4  := hfac$
iip4  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmp41 := coeffn(Pmmp4,eps,1)$
Qmmp41 := coeffn(Qmmp4,eps,1)$

Pmmp41_target := m * ((df(H4,r)-i*s*df(V4,r)) + (4+1)*(H4-i*s*V4)/r)/(2*4)$
Qmmp41_target := m * ((df(H4,r)-i*s*df(V4,r)) - (4-1)*(H4-i*s*V4)/r)/(2*4)$

write ("Pmmp41 residual:");
Pmmp41 - Pmmp41_target;
write ("Qmmp41 residual:");
Qmmp41 - Qmmp41_target;

P00p41 := coeffn(P00p4,eps,1)$

P00p41_target := (H4 - i*V4)/2/r$

write ("P00p41 residual:");
P00p41 - P00p41_target;

P01p41 := coeffn(P01p4,eps,1)$

P01p41_target := ((df(H4,r)-i*s*df(V4,r)) + (4+1)*(H4-i*s*V4)/r)/(2*4)$

write ("P01p41 residual:");
P01p41 - P01p41_target;

hhp41 := coeffn(hhp4,eps,1)$
iip41 := coeffn(iip4,eps,1)$

hhp41_target :=    df(H4,r) + i*s*df(V4,r)$
iip41_target := s*(df(H4,r) + i*s*df(V4,r))/4 + s*(4*4-1)*(H4 + i*s*V4)/4/r$

write ("hhp41 residual:");
hhp41 - hhp41_target;
write ("iip41 residual:");
iip41 - iip41_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=-1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = -1$
Pmmm1 := Pfac$
Qmmm1 := Qfac$
P00m1 := Pfac0$
P01m1 := Pfac1$
hhm1  := hfac$
iim1  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmm11 := coeffn(Pmmm1,eps,1)$
Qmmm11 := coeffn(Qmmm1,eps,1)$

Pmmm11_target := m * ( r/4/m - df(H1,r)/2)$
Qmmm11_target := m * ((1/4/m + 1/2)*r - df(H1,r)/2 - H1/r)$

write ("Pmmm11 residual:");
Pmmm11 - Pmmm11_target;
write ("Qmmm11 residual:");
Qmmm11 - Qmmm11_target;

P00m11 := coeffn(P00m1,eps,1)$

P00m11_target := lg*r/4 - r/4 + H1/2/r$

write ("P00m11 residual:");
P00m11 - P00m11_target;

P01m11 := coeffn(P01m1,eps,1)$

P01m11_target := (r/4 - df(H1,r)/2)$

write ("P01m11 residual:");
P01m11 - P01m11_target;

hhm11 := coeffn(hhm1,eps,1)$
iim11 := coeffn(iim1,eps,1)$

hhm11_target := - r + df(H1,r)$
iim11_target := - s*(-r + df(H1,r))$

write ("hhm11 residual:");
hhm11 - hhp11_target;
write ("iim11 residual:");
iim11 - iim11_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=-2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = -2$
Pmmm2 := Pfac$
Qmmm2 := Qfac$
P00m2 := Pfac0$
P01m2 := Pfac1$
hhm2  := hfac$
iim2  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmm21 := coeffn(Pmmm2,eps,1)$
Qmmm21 := coeffn(Qmmm2,eps,1)$

Pmmm21_target := m * ((- df(H2,r) - i*s*df(V2,r)) - (2-1)*(- H2 - i*s*V2)/r)/(2*2)$
Qmmm21_target := m * ((- df(H2,r) - i*s*df(V2,r)) + (2+1)*(- H2 - i*s*V2)/r)/(2*2)$

write ("Pmmm21 residual:");
Pmmm21 - Pmmm21_target;
write ("Qmmm21 residual:");
Qmmm21 - Qmmm21_target;

P00m21 := coeffn(P00m2,eps,1)$

P00m21_target := (H2 + i*V2)/2/r$

write ("P00m21 residual:");
P00m21 - P00m21_target;

P01m21 := coeffn(P01m2,eps,1)$

P01m21_target := ((- df(H2,r) - i*s*df(V2,r)) - (2-1)*(- H2 - i*s*V2)/r)/(2*2)$

write ("P01m21 residual:");
P01m21 - P01m21_target;

hhm21 := coeffn(hhm2,eps,1)$
iim21 := coeffn(iim2,eps,1)$

hhm21_target :=      df(H2,r) - i*s*df(V2,r)$
iim21_target := - s*(df(H2,r) - i*s*df(V2,r))/2 - s*(2*2-1)*(H2 - i*s*V2)/2/r$

write ("hhm21 residual:");
hhm21 - hhm21_target;
write ("iim21 residual:");
iim21 - iim21_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=-3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = -3$
Pmmm3 := Pfac$
Qmmm3 := Qfac$
P00m3 := Pfac0$
P01m3 := Pfac1$
hhm3  := hfac$
iim3  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmm31 := coeffn(Pmmm3,eps,1)$
Qmmm31 := coeffn(Qmmm3,eps,1)$

Pmmm31_target := m * ((- df(H3,r) - i*s*df(V3,r)) - (3-1)*(- H3 - i*s*V3)/r)/(2*3)$
Qmmm31_target := m * ((- df(H3,r) - i*s*df(V3,r)) + (3+1)*(- H3 - i*s*V3)/r)/(2*3)$

write ("Pmmm31 residual:");
Pmmm31 - Pmmm31_target;
write ("Qmmm31 residual:");
Qmmm31 - Qmmm31_target;

P00m31 := coeffn(P00m3,eps,1)$

P00m31_target := (H3 + i*V3)/2/r$

write ("P00m31 residual:");
P00m31 - P00m31_target;

P01m31 := coeffn(P01m3,eps,1)$

P01m31_target := ((- df(H3,r) - i*s*df(V3,r)) - (3-1)*(- H3 - i*s*V3)/r)/(2*3)$

write ("P01m31 residual:");
P01m31 - P01m31_target;

hhm31 := coeffn(hhm3,eps,1)$
iim31 := coeffn(iim3,eps,1)$

hhm31_target :=      df(H3,r) - i*s*df(V3,r)$
iim31_target := - s*(df(H3,r) - i*s*df(V3,r))/3 - s*(3*3-1)*(H3 - i*s*V3)/3/r$

write ("hhm31 residual:");
hhm31 - hhm31_target;
write ("iim31 residual:");
iim31 - iim31_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k=-4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let kk = -4$
Pmmm4 := Pfac$
Qmmm4 := Qfac$
P00m4 := Pfac0$
P01m4 := Pfac1$
hhm4  := hfac$
iim4  := ifac$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

Pmmm41 := coeffn(Pmmm4,eps,1)$
Qmmm41 := coeffn(Qmmm4,eps,1)$

Pmmm41_target := m * ((- df(H4,r) - i*s*df(V4,r)) - (4-1)*(- H4 - i*s*V4)/r)/(2*4)$
Qmmm41_target := m * ((- df(H4,r) - i*s*df(V4,r)) + (4+1)*(- H4 - i*s*V4)/r)/(2*4)$

write ("Pmmm41 residual:");
Pmmm41 - Pmmm41_target;
write ("Qmmm41 residual:");
Qmmm41 - Qmmm41_target;

P00m41 := coeffn(P00m4,eps,1)$

P00m41_target := (H4 + i*V4)/2/r$

write ("P00m41 residual:");
P00m41 - P00m41_target;

P01m41 := coeffn(P01m4,eps,1)$

P01m41_target :=((- df(H4,r) - i*s*df(V4,r)) - (4-1)*(- H4 - i*s*V4)/r)/(2*4)$

write ("P01m41 residual:");
P01m41 - P01m41_target;

hhm41 := coeffn(hhm4,eps,1)$
iim41 := coeffn(iim4,eps,1)$

hhm41_target :=      df(H4,r) - i*s*df(V4,r)$
iim41_target := - s*(df(H4,r) - i*s*df(V4,r))/4 - s*(4*4-1)*(H4 - i*s*V4)/4/r$

write ("hhmp4 residual:");
hhm41 - hhm41_target;
write ("iim41 residual:");
iim41 - iim41_target;

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate psi coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let df(H1,r,2) = H1pp$
let df(H2,r,2) = H2pp$
let df(H3,r,2) = H3pp$
let df(H4,r,2) = H4pp$
let df(V2,r,2) = V2pp$
let df(V3,r,2) = V3pp$
let df(V4,r,2) = V4pp$

Rmm0  := hh00_target *r**m  *r*df(Pmm0_target * r**(-m),r)$
Rmm2  :=
   + hh02_target  *r**m *r*df(Pmm0_target   * r**(-m),r)
   + hh00_target  *r**m *r*df(Pmm2_target   * r**(-m),r)
   + hhp11_target *r**m *r*df(Pmmp11_target * r**(-m),r)
   + hhp21_target *r**m *r*df(Pmmp21_target * r**(-m),r)
   + hhp31_target *r**m *r*df(Pmmp31_target * r**(-m),r)
   + hhp41_target *r**m *r*df(Pmmp41_target * r**(-m),r)
   + hhm11_target *r**m *r*df(Pmmm11_target * r**(-m),r)
   + hhm21_target *r**m *r*df(Pmmm21_target * r**(-m),r)
   + hhm31_target *r**m *r*df(Pmmm31_target * r**(-m),r)
   + hhm41_target *r**m *r*df(Pmmm41_target * r**(-m),r)
   + s*iip11_target *(m+1)*Pmmp11_target
   + s*iip21_target *(m+2)*Pmmp21_target
   + s*iip31_target *(m+3)*Pmmp31_target
   + s*iip41_target *(m+4)*Pmmp41_target
   + s*iim11_target *(m-1)*Pmmm11_target
   + s*iim21_target *(m-2)*Pmmm21_target
   + s*iim31_target *(m-3)*Pmmm31_target
   + s*iim41_target *(m-4)*Pmmm41_target$

G3 := - (((m - 2)*n**2 + m**2/4) /4/(m-1)/m)*r**2$

Rmm0_target := - m$
Rmm2_target := - m * (G3 - m*(G1 + H1/2) - m**2*G2)$

write ("Rmm0 residual:");
Rmm0 - Rmm0_target;

write ("Rmm2 residual:");
Rmm2 - Rmm2_target;

R000 := hh00_target *r*df(P000_target,r)$

R002 := 
   + hh02_target  *r*df(P000_target  ,r)
   + hh00_target  *r*df(P002_target  ,r)
   + hhp11_target *r*df(P00p11_target,r)
   + hhp21_target *r*df(P00p21_target,r)
   + hhp31_target *r*df(P00p31_target,r)
   + hhp41_target *r*df(P00p41_target,r)
   + hhm11_target *r*df(P00m11_target,r)
   + hhm21_target *r*df(P00m21_target,r)
   + hhm31_target *r*df(P00m31_target,r)
   + hhm41_target *r*df(P00m41_target,r)
   + iip11_target *1   *P00p11_target
   + iip21_target *2   *P00p21_target
   + iip31_target *3   *P00p31_target
   + iip41_target *4   *P00p41_target
   + iim11_target *(-1)*P00m11_target
   + iim21_target *(-2)*P00m21_target
   + iim31_target *(-3)*P00m31_target
   + iim41_target *(-4)*P00m41_target$
let df(lg,r) = -1/r$
R002 := sub(s=1,R002)$

G3b := - n**2*(lg/2 + 1/4)*r**2$

R000_target := - 1$
R002_target := - G3b$

write ("R000 residual:");
R000 - R000_target;

write ("R002 residual:");
R002 - R002_target;

R010 := hh00_target *r *r*df(P010_target/r,r)$

R012 := 
   + hh02_target  *r *r*df(P010_target  /r,r)
   + hh00_target  *r *r*df(P012_target  /r,r)
   + hhp11_target *r *r*df(P01p11_target/r,r)
   + hhp21_target *r *r*df(P01p21_target/r,r)
   + hhp31_target *r *r*df(P01p31_target/r,r)
   + hhp41_target *r *r*df(P01p41_target/r,r)
   + hhm11_target *r *r*df(P01m11_target/r,r)
   + hhm21_target *r *r*df(P01m21_target/r,r)
   + hhm31_target *r *r*df(P01m31_target/r,r)
   + hhm41_target *r *r*df(P01m41_target/r,r)
   + s*iip11_target *(1+1)*P01p11_target
   + s*iip21_target *(1+2)*P01p21_target
   + s*iip31_target *(1+3)*P01p31_target
   + s*iip41_target *(1+4)*P01p41_target
   + s*iim11_target *(1-1)*P01m11_target
   + s*iim21_target *(1-2)*P01m21_target
   + s*iim31_target *(1-3)*P01m31_target
   + s*iim41_target *(1-4)*P01m41_target$

G3c := - (n**2/4 - 1/8 - (n**2/2 - 1/8)*lg)*r**2$

R010_target := - 1$
R012_target := - (G3c - (G1 + H1/2) - G2)$

write ("R010 residual:");
R010 - R010_target;

write ("R012 residual:");
R012 - R012_target;

Rmmp1 :=
     hh00_target    *r**m *r*df(Pmmp11_target * r**(-m),r)
   + hhm11_target   *r**m *r*df(Pmm0_target   * r**(-m),r)
   + s*iim11_target *m*Pmm0_target$

Rmmp1_target := - m * (1+m) * (- (1/4/m + 1/2)*r + H1/r + df(H1,r)/2)$

write ("Rmmp1 residual:");
Rmmp1 - Rmmp1_target;

R00p1 :=
     hh00_target  *r*df(P00p11_target,r)
   + hhm11_target *r*df(P000_target  ,r)$

R00p1_target := r*lg/4 + r/2 - H1/2/r - df(H1,r)/2$

write ("R00p1 residual:");
R00p1 - R00p1_target;

Rmmp2 :=
     hh00_target    *r**m *r*df(Pmmp21_target *r**(-m),r)
   + hhm21_target   *r**m *r*df(Pmm0_target   *r**(-m),r)
   + s*iim21_target *m*Pmm0_target$

Rmmp2_target := - m * (1+m/2) * ((df(H2,r) - i*s*df(V2,r))/2 + ((2+1)/2)*(H2 - i*s*V2)/r)$

write ("Rmmp2 residual:");
Rmmp2 - Rmmp2_target;

R00p2 :=
       hh00_target  *r*df(P00p21_target,r)
     + hhm21_target *r*df(P000_target  ,r)$
R00p2 := sub(s=1,R00p2)$

R00p2_target := - (df(H2,r) - i*df(V2,r))/2 - (H2 - i*V2)/2/r$

write ("Rmmp2 residual:");
Rmmp2 - Rmmp2_target;

Rmmp3 :=
     hh00_target    *r**m *r*df(Pmmp31_target *r**(-m),r)
   + hhm31_target   *r**m *r*df(Pmm0_target   *r**(-m),r)
   + s*iim31_target *m*Pmm0_target$

Rmmp3_target := - m * (1+m/3) * ((df(H3,r) - i*s*df(V3,r))/2 + ((3+1)/2)*(H3 - i*s*V3)/r)$

write ("Rmmp3 residual:");
Rmmp3 - Rmmp3_target;

R00p3 :=
       hh00_target  *r*df(P00p31_target,r)
     + hhm31_target *r*df(P000_target  ,r)$
R00p3 := sub(s=1,R00p3)$

R00p3_target := - (df(H3,r) - i*df(V3,r))/2 - (H3 - i*V3)/2/r$

write ("Rmmp3 residual:");
Rmmp3 - Rmmp3_target;

Rmmp4 :=
     hh00_target    *r**m *r*df(Pmmp41_target *r**(-m),r)
   + hhm41_target   *r**m *r*df(Pmm0_target   *r**(-m),r)
   + s*iim41_target *m*Pmm0_target$

Rmmp4_target := - m * (1+m/4) * ((df(H4,r) - i*s*df(V4,r))/2 + ((4+1)/2)*(H4 - i*s*V4)/r)$

write ("Rmmp4 residual:");
Rmmp4 - Rmmp4_target;

R00p4 :=
       hh00_target  *r*df(P00p41_target,r)
     + hhm41_target *r*df(P000_target  ,r)$
R00p4 := sub(s=1,R00p4)$

R00p4_target := - (df(H4,r) - i*df(V4,r))/2 - (H4 - i*V4)/2/r$

write ("Rmmp4 residual:");
Rmmp4 - Rmmp4_target;

Rmmm1 :=
     hh00_target    *r**m *r*df(Pmmm11_target * r**(-m),r)
   + hhp11_target   *r**m *r*df(Pmm0_target   * r**(-m),r)
   + s*iip11_target *m*Pmm0_target$

Rmmm1_target := - m * (- (1/4/m + 1/4)*r + (1 - m) * df(H1,r)/2)$

write ("Rmmm1 residual:");
Rmmm1 - Rmmm1_target;

R00m1 :=
     hh00_target  *r*df(P00m11_target,r)
   + hhp11_target *r*df(P000_target  ,r)$

R00m1_target := r*lg/4 + r/2 - H1/2/r - df(H1,r)/2$

write ("R00m1 residual:");
R00m1 - R00m1_target;

Rmmm2 :=
     hh00_target    *r**m *r*df(Pmmm21_target *r**(-m),r)
   + hhp21_target   *r**m *r*df(Pmm0_target   *r**(-m),r)
   + s*iip21_target *m*Pmm0_target$

Rmmm2_target := - m * (1-m/2) * ((df(H2,r) + i*s*df(V2,r))/2 - ((2-1)/2)*(H2 + i*s*V2)/r)$

write ("Rmmm2 residual:");
Rmmm2 - Rmmm2_target;

R00m2 :=
       hh00_target  *r*df(P00m21_target,r)
     + hhp21_target *r*df(P000_target  ,r)$
R00m2 := sub(s=1,R00m2)$

R00m2_target := - (df(H2,r) + i*df(V2,r))/2 - (H2 + i*V2)/2/r$

write ("R00m2 residual:");
R00m2 - R00m2_target;

Rmmm3 :=
     hh00_target    *r**m *r*df(Pmmm31_target *r**(-m),r)
   + hhp31_target   *r**m *r*df(Pmm0_target   *r**(-m),r)
   + s*iip31_target *m*Pmm0_target$

Rmmm3_target := - m * (1-m/3) * ((df(H3,r) + i*s*df(V3,r))/2 - ((3-1)/2)*(H3 + i*s*V3)/r)$

write ("Rmmm3 residual:");
Rmmm3 - Rmmm3_target;

R00m3 :=
       hh00_target  *r*df(P00m31_target,r)
     + hhp31_target *r*df(P000_target  ,r)$
R00m3 := sub(s=1,R00m3)$

R00m3_target := - (df(H3,r) + i*df(V3,r))/2 - (H3 + i*V3)/2/r$

write ("R00m3 residual:");
R00m3 - R00m3_target;

Rmmm4 :=
     hh00_target    *r**m *r*df(Pmmm41_target *r**(-m),r)
   + hhp41_target   *r**m *r*df(Pmm0_target   *r**(-m),r)
   + s*iip41_target *m*Pmm0_target$

Rmmm4_target := - m * (1-m/4) * ((df(H4,r) + i*s*df(V4,r))/2 - ((4-1)/2)*(H4 + i*s*V4)/r)$

write ("Rmmm4 residual:");
Rmmm4 - Rmmm4_target;

R00m4 :=
       hh00_target  *r*df(P00m41_target,r)
     + hhp41_target *r*df(P000_target  ,r)$
R00m2 := sub(s=1,R00m4)$

R00m4_target := - (df(H4,r) + i*df(V4,r))/2 - (H4 + i*V4)/2/r$

write ("R00m3 residual:");
R00m3 - R00m3_target;

Smm0  := hh00_target *r**(-m)  *r*df(Qmm0_target * r**(m),r)$
Smm2  :=
   + hh02_target  *r**(-m) *r*df(Qmm0_target   * r**(m),r)
   + hh00_target  *r**(-m) *r*df(Qmm2_target   * r**(m),r)
   + hhp11_target *r**(-m) *r*df(Qmmp11_target * r**(m),r)
   + hhp21_target *r**(-m) *r*df(Qmmp21_target * r**(m),r)
   + hhp31_target *r**(-m) *r*df(Qmmp31_target * r**(m),r)
   + hhp41_target *r**(-m) *r*df(Qmmp41_target * r**(m),r)
   + hhm11_target *r**(-m) *r*df(Qmmm11_target * r**(m),r)
   + hhm21_target *r**(-m) *r*df(Qmmm21_target * r**(m),r)
   + hhm31_target *r**(-m) *r*df(Qmmm31_target * r**(m),r)
   + hhm41_target *r**(-m) *r*df(Qmmm41_target * r**(m),r)
   + s*iip11_target *(m+1)*Qmmp11_target
   + s*iip21_target *(m+2)*Qmmp21_target
   + s*iip31_target *(m+3)*Qmmp31_target
   + s*iip41_target *(m+4)*Qmmp41_target
   + s*iim11_target *(m-1)*Qmmm11_target
   + s*iim21_target *(m-2)*Qmmm21_target
   + s*iim31_target *(m-3)*Qmmm31_target
   + s*iim41_target *(m-4)*Qmmm41_target$

G3a := (((m + 2)*n**2 - m**2/4) /4/(m+1)/m)*r**2$

Smm0_target := m$
Smm2_target := m * (G3a + m*(G1 + H1/2) - m**2*G2a)$

write ("Smm0 residual:");
Smm0 - Smm0_target;

write ("Smm2 residual:");
Smm2 - Smm2_target;

Smmp1 :=
     hh00_target    *r**(-m) *r*df(Qmmp11_target * r**(m),r)
   + hhm11_target   *r**(-m) *r*df(Qmm0_target   * r**(m),r)
   + s*iim11_target *m*Qmm0_target$

Smmp1_target := m * (- (-1/4/m + 1/4)*r + (1 + m) * df(H1,r)/2)$

write ("Smmp1 residual:");
Smmp1 - Smmp1_target;

Smmp2 :=
     hh00_target    *r**(-m) *r*df(Qmmp21_target *r**(m),r)
   + hhm21_target   *r**(-m) *r*df(Qmm0_target   *r**(m),r)
   + s*iim21_target *m*Qmm0_target$

Smmp2_target := m * (1+m/2) * ((df(H2,r) - i*s*df(V2,r))/2 - ((2-1)/2)*(H2 - i*s*V2)/r)$

write ("Smmp2 residual:");
Smmp2 - Smmp2_target;

Smmp3 :=
     hh00_target    *r**(-m) *r*df(Qmmp31_target *r**(m),r)
   + hhm31_target   *r**(-m) *r*df(Qmm0_target   *r**(m),r)
   + s*iim31_target *m*Qmm0_target$

Smmp3_target := m * (1+m/3) * ((df(H3,r) - i*s*df(V3,r))/2 - ((3-1)/2)*(H3 - i*s*V3)/r)$

write ("Smmp3 residual:");
Smmp3 - Smmp3_target;

Smmp4 :=
     hh00_target    *r**(-m) *r*df(Qmmp41_target *r**(m),r)
   + hhm41_target   *r**(-m) *r*df(Qmm0_target   *r**(m),r)
   + s*iim41_target *m*Qmm0_target$

Smmp4_target := m * (1+m/4) * ((df(H4,r) - i*s*df(V4,r))/2 - ((4-1)/2)*(H4 - i*s*V4)/r)$

write ("Smmp4 residual:");
Smmp4 - Smmp4_target;

Smmm1 :=
     hh00_target    *r**(-m) *r*df(Qmmm11_target * r**(m),r)
   + hhp11_target   *r**(-m) *r*df(Qmm0_target   * r**(m),r)
   + s*iip11_target *m*Qmm0_target$

Smmm1_target := m * (1-m) * (- (-1/4/m + 1/2)*r + H1/r + df(H1,r)/2)$

write ("Smmm1 residual:");
Smmm1 - Smmm1_target;

Smmm2 :=
     hh00_target    *r**(-m) *r*df(Qmmm21_target * r**(m),r)
   + hhp21_target   *r**(-m) *r*df(Qmm0_target   * r**(m),r)
   + s*iip21_target *m*Qmm0_target$

Smmm2_target := m * (1-m/2) * ((df(H2,r) + i*s*df(V2,r))/2 + ((2+1)/2)*(H2 +i*s*V2)/r)$

write ("Smmm2 residual:");
Smmm2 - Smmm2_target;

Smmm3 :=
     hh00_target    *r**(-m) *r*df(Qmmm31_target * r**(m),r)
   + hhp31_target   *r**(-m) *r*df(Qmm0_target   * r**(m),r)
   + s*iip31_target *m*Qmm0_target$

Smmm3_target := m * (1-m/3) * ((df(H3,r) + i*s*df(V3,r))/2 + ((3+1)/2)*(H3 +i*s*V3)/r)$

write ("Smmm3 residual:");
Smmm3 - Smmm3_target;

Smmm4 :=
     hh00_target    *r**(-m) *r*df(Qmmm41_target * r**(m),r)
   + hhp41_target   *r**(-m) *r*df(Qmm0_target   * r**(m),r)
   + s*iip41_target *m*Qmm0_target$

Smmm4_target := m * (1-m/4) * ((df(H4,r) + i*s*df(V4,r))/2 + ((4+1)/2)*(H4 +i*s*V4)/r)$

write ("Smmm4 residual:");
Smmm4 - Smmm4_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to change mode number in expression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure sft(ex,k)$
   begin scalar c$
      c := sub(m=j,ex)$
      c := sub(j=m+k,c)$
      return c$
   end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate general elements of solution vectors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pmm := Pmm0_target + eps**2*Pmm2_target$
Qmm := Qmm0_target + eps**2*Qmm2_target$
Rmm := Rmm0_target + eps**2*Rmm2_target$
Smm := Smm0_target + eps**2*Smm2_target$

Pmp1m := eps * Pmmp11_target$
Qmp1m := eps * Qmmp11_target$
Rmp1m := eps * Rmmp1_target$
Smp1m := eps * Smmp1_target$

Pmp2m := eps * Pmmp21_target$
Qmp2m := eps * Qmmp21_target$
Rmp2m := eps * Rmmp2_target$
Smp2m := eps * Smmp2_target$

Pmp3m := eps * Pmmp31_target$
Qmp3m := eps * Qmmp31_target$
Rmp3m := eps * Rmmp3_target$
Smp3m := eps * Smmp3_target$

Pmp4m := eps * Pmmp41_target$
Qmp4m := eps * Qmmp41_target$
Rmp4m := eps * Rmmp4_target$
Smp4m := eps * Smmp4_target$

Pmm1m := eps * Pmmm11_target$
Qmm1m := eps * Qmmm11_target$
Rmm1m := eps * Rmmm1_target$
Smm1m := eps * Smmm1_target$

Pmm2m := eps * Pmmm21_target$
Qmm2m := eps * Qmmm21_target$
Rmm2m := eps * Rmmm2_target$
Smm2m := eps * Smmm2_target$

Pmm3m := eps * Pmmm31_target$
Qmm3m := eps * Qmmm31_target$
Rmm3m := eps * Rmmm3_target$
Smm3m := eps * Smmm3_target$

Pmm4m := eps * Pmmm41_target$
Qmm4m := eps * Qmmm41_target$
Rmm4m := eps * Rmmm4_target$
Smm4m := eps * Smmm4_target$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amm :=
   + sub(i=-i,Pmm)  *Rmm   - sub(i=-i,Rmm)  *Pmm
   + sub(i=-i,Pmp1m)*Rmp1m - sub(i=-i,Rmp1m)*Pmp1m
   + sub(i=-i,Pmp2m)*Rmp2m - sub(i=-i,Rmp2m)*Pmp2m
   + sub(i=-i,Pmp3m)*Rmp3m - sub(i=-i,Rmp3m)*Pmp3m
   + sub(i=-i,Pmp4m)*Rmp4m - sub(i=-i,Rmp4m)*Pmp4m
   + sub(i=-i,Pmm1m)*Rmm1m - sub(i=-i,Rmm1m)*Pmm1m
   + sub(i=-i,Pmm2m)*Rmm2m - sub(i=-i,Rmm2m)*Pmm2m
   + sub(i=-i,Pmm3m)*Rmm3m - sub(i=-i,Rmm3m)*Pmm3m
   + sub(i=-i,Pmm4m)*Rmm4m - sub(i=-i,Rmm4m)*Pmm4m;

Bmm :=
   + sub(i=-i,Pmm)  *Smm   - sub(i=-i,Rmm)  *Qmm
   + sub(i=-i,Pmp1m)*Smp1m - sub(i=-i,Rmp1m)*Qmp1m
   + sub(i=-i,Pmp2m)*Smp2m - sub(i=-i,Rmp2m)*Qmp2m
   + sub(i=-i,Pmp3m)*Smp3m - sub(i=-i,Rmp3m)*Qmp3m
   + sub(i=-i,Pmp4m)*Smp4m - sub(i=-i,Rmp4m)*Qmp4m
   + sub(i=-i,Pmm1m)*Smm1m - sub(i=-i,Rmm1m)*Qmm1m
   + sub(i=-i,Pmm2m)*Smm2m - sub(i=-i,Rmm2m)*Qmm2m
   + sub(i=-i,Pmm3m)*Smm3m - sub(i=-i,Rmm3m)*Qmm3m
   + sub(i=-i,Pmm4m)*Smm4m - sub(i=-i,Rmm4m)*Qmm4m;

Cmm :=
   + sub(i=-i,Qmm)  *Smm   - sub(i=-i,Smm)  *Qmm
   + sub(i=-i,Qmp1m)*Smp1m - sub(i=-i,Smp1m)*Qmp1m
   + sub(i=-i,Qmp2m)*Smp2m - sub(i=-i,Smp2m)*Qmp2m
   + sub(i=-i,Qmp3m)*Smp3m - sub(i=-i,Smp3m)*Qmp3m
   + sub(i=-i,Qmp4m)*Smp4m - sub(i=-i,Smp4m)*Qmp4m
   + sub(i=-i,Qmm1m)*Smm1m - sub(i=-i,Smm1m)*Qmm1m
   + sub(i=-i,Qmm2m)*Smm2m - sub(i=-i,Smm2m)*Qmm2m
   + sub(i=-i,Qmm3m)*Smm3m - sub(i=-i,Smm3m)*Qmm3m
   + sub(i=-i,Qmm4m)*Smm4m - sub(i=-i,Smm4m)*Qmm4m;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m+1 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ammp1 :=
   + sub(i=-i,Pmm)  *sft(Rmm1m,1) - sub(i=-i,Rmm)  *sft(Pmm1m,1)
   + sub(i=-i,Pmp1m)*sft(Rmm,  1) - sub(i=-i,Rmp1m)*sft(Pmm,  1)
   + sub(i=-i,Pmp2m)*sft(Rmp1m,1) - sub(i=-i,Rmp2m)*sft(Pmp2m,1)
   + sub(i=-i,Pmp3m)*sft(Rmp2m,1) - sub(i=-i,Rmp3m)*sft(Pmp3m,1)
   + sub(i=-i,Pmp4m)*sft(Rmp3m,1) - sub(i=-i,Rmp4m)*sft(Pmp4m,1)
   + sub(i=-i,Pmm1m)*sft(Rmm2m,1) - sub(i=-i,Rmm1m)*sft(Pmm2m,1)
   + sub(i=-i,Pmm2m)*sft(Rmm3m,1) - sub(i=-i,Rmm2m)*sft(Pmm3m,1)
   + sub(i=-i,Pmm3m)*sft(Rmm4m,1) - sub(i=-i,Rmm3m)*sft(Pmm4m,1)$
Ammp1 := coeffn(Ammp1,eps,0) + eps*coeffn(Ammp1,eps,1);

Bmmp1 :=
   + sub(i=-i,Pmm)  *sft(Smm1m,1) - sub(i=-i,Rmm)  *sft(Qmm1m,1)
   + sub(i=-i,Pmp1m)*sft(Smm,  1) - sub(i=-i,Rmp1m)*sft(Qmm,  1)
   + sub(i=-i,Pmp2m)*sft(Smp1m,1) - sub(i=-i,Rmp2m)*sft(Qmp2m,1)
   + sub(i=-i,Pmp3m)*sft(Smp2m,1) - sub(i=-i,Rmp3m)*sft(Qmp3m,1)
   + sub(i=-i,Pmp4m)*sft(Smp3m,1) - sub(i=-i,Rmp4m)*sft(Qmp4m,1)
   + sub(i=-i,Pmm1m)*sft(Smm2m,1) - sub(i=-i,Rmm1m)*sft(Qmm2m,1)
   + sub(i=-i,Pmm2m)*sft(Smm3m,1) - sub(i=-i,Rmm2m)*sft(Qmm3m,1)
   + sub(i=-i,Pmm3m)*sft(Smm4m,1) - sub(i=-i,Rmm3m)*sft(Qmm4m,1)$
Bmmp1 := coeffn(Bmmp1,eps,0) + eps*coeffn(Bmmp1,eps,1);

Cmmp1 :=
   + sub(i=-i,Qmm)  *sft(Smm1m,1) - sub(i=-i,Smm)  *sft(Qmm1m,1)
   + sub(i=-i,Qmp1m)*sft(Smm,  1) - sub(i=-i,Smp1m)*sft(Qmm,  1)
   + sub(i=-i,Qmp2m)*sft(Smp1m,1) - sub(i=-i,Smp2m)*sft(Qmp2m,1)
   + sub(i=-i,Qmp3m)*sft(Smp2m,1) - sub(i=-i,Smp3m)*sft(Qmp3m,1)
   + sub(i=-i,Qmp4m)*sft(Smp3m,1) - sub(i=-i,Smp4m)*sft(Qmp4m,1)
   + sub(i=-i,Qmm1m)*sft(Smm2m,1) - sub(i=-i,Smm1m)*sft(Qmm2m,1)
   + sub(i=-i,Qmm2m)*sft(Smm3m,1) - sub(i=-i,Smm2m)*sft(Qmm3m,1)
   + sub(i=-i,Qmm3m)*sft(Smm4m,1) - sub(i=-i,Smm3m)*sft(Qmm4m,1)$
Cmmp1 := coeffn(Cmmp1,eps,0) + eps*coeffn(Cmmp1,eps,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m+2 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ammp2 :=
   + sub(i=-i,Pmm)  *sft(Rmm2m,2) - sub(i=-i,Rmm)  *sft(Pmm2m,2)
   + sub(i=-i,Pmp1m)*sft(Rmm1m,2) - sub(i=-i,Rmp1m)*sft(Pmm1m,2)
   + sub(i=-i,Pmp2m)*sft(Rmm,  2) - sub(i=-i,Rmp2m)*sft(Pmm,  2)
   + sub(i=-i,Pmp3m)*sft(Rmp1m,2) - sub(i=-i,Rmp3m)*sft(Pmp1m,2)
   + sub(i=-i,Pmp4m)*sft(Rmp2m,2) - sub(i=-i,Rmp4m)*sft(Pmp2m,2)
   + sub(i=-i,Pmm1m)*sft(Rmm3m,2) - sub(i=-i,Rmm1m)*sft(Pmm3m,2)
   + sub(i=-i,Pmm2m)*sft(Rmm4m,2) - sub(i=-i,Rmm2m)*sft(Pmm4m,2)$
Ammp2 := coeffn(Ammp2,eps,0) + eps*coeffn(Ammp2,eps,1);

Bmmp2 :=
   + sub(i=-i,Pmm)  *sft(Smm2m,2) - sub(i=-i,Rmm)  *sft(Qmm2m,2)
   + sub(i=-i,Pmp1m)*sft(Smm1m,2) - sub(i=-i,Rmp1m)*sft(Qmm1m,2)
   + sub(i=-i,Pmp2m)*sft(Smm,  2) - sub(i=-i,Rmp2m)*sft(Qmm,  2)
   + sub(i=-i,Pmp3m)*sft(Smp1m,2) - sub(i=-i,Rmp3m)*sft(Qmp1m,2)
   + sub(i=-i,Pmp4m)*sft(Smp2m,2) - sub(i=-i,Rmp4m)*sft(Qmp2m,2)
   + sub(i=-i,Pmm1m)*sft(Smm3m,2) - sub(i=-i,Rmm1m)*sft(Qmm3m,2)
   + sub(i=-i,Pmm2m)*sft(Smm4m,2) - sub(i=-i,Rmm2m)*sft(Qmm4m,2)$
Bmmp2 := coeffn(Bmmp2,eps,0) + eps*coeffn(Bmmp2,eps,1);

Cmmp2 :=
   + sub(i=-i,Qmm)  *sft(Smm2m,2) - sub(i=-i,Smm)  *sft(Qmm2m,2)
   + sub(i=-i,Qmp1m)*sft(Smm1m,2) - sub(i=-i,Smp1m)*sft(Qmm1m,2)
   + sub(i=-i,Qmp2m)*sft(Smm,  2) - sub(i=-i,Smp2m)*sft(Qmm,  2)
   + sub(i=-i,Qmp3m)*sft(Smp1m,2) - sub(i=-i,Smp3m)*sft(Qmp1m,2)
   + sub(i=-i,Qmp4m)*sft(Smp2m,2) - sub(i=-i,Smp4m)*sft(Qmp2m,2)
   + sub(i=-i,Qmm1m)*sft(Smm3m,2) - sub(i=-i,Smm1m)*sft(Qmm3m,2)
   + sub(i=-i,Qmm2m)*sft(Smm4m,2) - sub(i=-i,Smm2m)*sft(Qmm4m,2)$
 Cmmp2 := coeffn(Cmmp2,eps,0) + eps*coeffn(Cmmp2,eps,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m+3 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ammp3 :=
   + sub(i=-i,Pmm)  *sft(Rmm3m,3) - sub(i=-i,Rmm)  *sft(Pmm3m,3)
   + sub(i=-i,Pmp1m)*sft(Rmm2m,3) - sub(i=-i,Rmp1m)*sft(Pmm2m,3)
   + sub(i=-i,Pmp2m)*sft(Rmm1m,3) - sub(i=-i,Rmp2m)*sft(Pmm1m,3)
   + sub(i=-i,Pmp3m)*sft(Rmm,  3) - sub(i=-i,Rmp3m)*sft(Pmm,  3)
   + sub(i=-i,Pmp4m)*sft(Rmp1m,3) - sub(i=-i,Rmp4m)*sft(Pmp1m,3)
   + sub(i=-i,Pmm1m)*sft(Rmm4m,3) - sub(i=-i,Rmm1m)*sft(Pmm4m,3)$
Ammp3 := coeffn(Ammp3,eps,0) + eps*coeffn(Ammp3,eps,1);

Bmmp3 :=
   + sub(i=-i,Pmm)  *sft(Smm3m,3) - sub(i=-i,Rmm)  *sft(Qmm3m,3)
   + sub(i=-i,Pmp1m)*sft(Smm2m,3) - sub(i=-i,Rmp1m)*sft(Qmm2m,3)
   + sub(i=-i,Pmp2m)*sft(Smm1m,3) - sub(i=-i,Rmp2m)*sft(Qmm1m,3)
   + sub(i=-i,Pmp3m)*sft(Smm,  3) - sub(i=-i,Rmp3m)*sft(Qmm,  3)
   + sub(i=-i,Pmp4m)*sft(Rmp1m,3) - sub(i=-i,Rmp4m)*sft(Qmp1m,3)
   + sub(i=-i,Pmm1m)*sft(Smm4m,3) - sub(i=-i,Rmm1m)*sft(Qmm4m,3)$
Bmmp3 := coeffn(Bmmp3,eps,0) + eps*coeffn(Bmmp3,eps,1);

Cmmp3 :=
   + sub(i=-i,Qmm)  *sft(Smm3m,3) - sub(i=-i,Smm)  *sft(Qmm3m,3)
   + sub(i=-i,Qmp1m)*sft(Smm3m,3) - sub(i=-i,Smp1m)*sft(Qmm2m,3)
   + sub(i=-i,Qmp2m)*sft(Smm1m,3) - sub(i=-i,Smp2m)*sft(Qmm1m,3)
   + sub(i=-i,Qmp3m)*sft(Smm,  3) - sub(i=-i,Smp3m)*sft(Qmm,  3)
   + sub(i=-i,Qmp4m)*sft(Smp1m,3) - sub(i=-i,Smp4m)*sft(Qmp1m,3)
   + sub(i=-i,Qmm1m)*sft(Smm4m,3) - sub(i=-i,Smm1m)*sft(Qmm2m,3)$
Cmmp3 := coeffn(Cmmp3,eps,0) + eps*coeffn(Cmmp3,eps,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m-1 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ammm1 :=
   + sub(i=-i,Pmm)  *sft(Rmp1m,-1) - sub(i=-i,Rmm)  *sft(Pmp1m,-1)
   + sub(i=-i,Pmp1m)*sft(Rmp2m,-1) - sub(i=-i,Rmp1m)*sft(Pmp2m,-1)
   + sub(i=-i,Pmp2m)*sft(Rmp3m,-1) - sub(i=-i,Rmp2m)*sft(Pmp3m,-1)
   + sub(i=-i,Pmp3m)*sft(Rmp4m,-1) - sub(i=-i,Rmp3m)*sft(Pmp4m,-1)
   + sub(i=-i,Pmm1m)*sft(Rmm,  -1) - sub(i=-i,Rmm1m)*sft(Pmm,  -1)
   + sub(i=-i,Pmm2m)*sft(Rmm1m,-1) - sub(i=-i,Rmm2m)*sft(Pmm1m,-1)
   + sub(i=-i,Pmm3m)*sft(Rmm2m,-1) - sub(i=-i,Rmm3m)*sft(Pmm2m,-1)
   + sub(i=-i,Pmm4m)*sft(Rmm3m,-1) - sub(i=-i,Rmm4m)*sft(Pmm3m,-1)$
Ammm1 := coeffn(Ammm1,eps,0) + eps*coeffn(Ammm1,eps,1);

Bmmm1 :=
   + sub(i=-i,Pmm)  *sft(Smp1m,-1) - sub(i=-i,Rmm)  *sft(Qmp1m,-1)
   + sub(i=-i,Pmp1m)*sft(Smp2m,-1) - sub(i=-i,Rmp1m)*sft(Qmp2m,-1)
   + sub(i=-i,Pmp2m)*sft(Smp3m,-1) - sub(i=-i,Rmp2m)*sft(Qmp3m,-1)
   + sub(i=-i,Pmp3m)*sft(Smp4m,-1) - sub(i=-i,Rmp3m)*sft(Qmp4m,-1)
   + sub(i=-i,Pmm1m)*sft(Smm,  -1) - sub(i=-i,Rmm1m)*sft(Qmm,  -1)
   + sub(i=-i,Pmm2m)*sft(Smm1m,-1) - sub(i=-i,Rmm2m)*sft(Qmm1m,-1)
   + sub(i=-i,Pmm3m)*sft(Smm2m,-1) - sub(i=-i,Rmm3m)*sft(Qmm2m,-1)
   + sub(i=-i,Pmm4m)*sft(Smm3m,-1) - sub(i=-i,Rmm4m)*sft(Qmm3m,-1)$
Bmmm1 := coeffn(Bmmm1,eps,0) + eps*coeffn(Bmmm1,eps,1);

Cmmm1 :=
   + sub(i=-i,Qmm)  *sft(Smp1m,-1) - sub(i=-i,Smm)  *sft(Qmp1m,-1)
   + sub(i=-i,Qmp1m)*sft(Smp2m,-1) - sub(i=-i,Smp1m)*sft(Qmp2m,-1)
   + sub(i=-i,Qmp2m)*sft(Smp3m,-1) - sub(i=-i,Smp2m)*sft(Qmp3m,-1)
   + sub(i=-i,Qmp3m)*sft(Smp4m,-1) - sub(i=-i,Smp3m)*sft(Qmp4m,-1)
   + sub(i=-i,Qmm1m)*sft(Smm,  -1) - sub(i=-i,Smm1m)*sft(Qmm,  -1)
   + sub(i=-i,Qmm2m)*sft(Smm1m,-1) - sub(i=-i,Smm2m)*sft(Qmm1m,-1)
   + sub(i=-i,Qmm3m)*sft(Smm2m,-1) - sub(i=-i,Smm3m)*sft(Qmm2m,-1)
   + sub(i=-i,Qmm4m)*sft(Smm3m,-1) - sub(i=-i,Smm4m)*sft(Qmm3m,-1)$
Cmmm1 := coeffn(Cmmm1,eps,0) + eps*coeffn(Cmmm1,eps,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m-2 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ammm2 :=
   + sub(i=-i,Pmm)  *sft(Rmp2m,-2) - sub(i=-i,Rmm)  *sft(Pmp2m,-2)
   + sub(i=-i,Pmp1m)*sft(Rmp3m,-2) - sub(i=-i,Rmp1m)*sft(Pmp3m,-2)
   + sub(i=-i,Pmp2m)*sft(Rmp4m,-2) - sub(i=-i,Rmp2m)*sft(Pmp4m,-2)
   + sub(i=-i,Pmm1m)*sft(Rmp1m,-2) - sub(i=-i,Rmm1m)*sft(Pmp1m,-2)
   + sub(i=-i,Pmm2m)*sft(Rmm,  -2) - sub(i=-i,Rmm2m)*sft(Pmm,  -2)
   + sub(i=-i,Pmm3m)*sft(Rmm1m,-2) - sub(i=-i,Rmm3m)*sft(Pmm1m,-2)
   + sub(i=-i,Pmm4m)*sft(Rmm2m,-2) - sub(i=-i,Rmm4m)*sft(Pmm2m,-2)$
Ammm2 := coeffn(Ammm2,eps,0) + eps*coeffn(Ammm2,eps,1);

Bmmm2 :=
   + sub(i=-i,Pmm)  *sft(Smp2m,-2) - sub(i=-i,Rmm)  *sft(Qmp2m,-2)
   + sub(i=-i,Pmp1m)*sft(Smp3m,-2) - sub(i=-i,Rmp1m)*sft(Qmp3m,-2)
   + sub(i=-i,Pmp2m)*sft(Smp4m,-2) - sub(i=-i,Rmp2m)*sft(Qmp4m,-2)
   + sub(i=-i,Pmm1m)*sft(Smp1m,-2) - sub(i=-i,Rmm1m)*sft(Qmp1m,-2)
   + sub(i=-i,Pmm2m)*sft(Smm,  -2) - sub(i=-i,Rmm2m)*sft(Qmm,  -2)
   + sub(i=-i,Pmm3m)*sft(Smm1m,-2) - sub(i=-i,Rmm3m)*sft(Qmm1m,-2)
   + sub(i=-i,Pmm4m)*sft(Smm2m,-2) - sub(i=-i,Rmm4m)*sft(Qmm2m,-2)$
Bmmm2 := coeffn(Bmmm2,eps,0) + eps*coeffn(Bmmm2,eps,1);

Cmmm2 :=
   + sub(i=-i,Qmm)  *sft(Smp2m,-2) - sub(i=-i,Smm)  *sft(Qmp2m,-2)
   + sub(i=-i,Qmp1m)*sft(Smp3m,-2) - sub(i=-i,Smp1m)*sft(Qmp3m,-2)
   + sub(i=-i,Qmp2m)*sft(Smp4m,-2) - sub(i=-i,Smp2m)*sft(Qmp4m,-2)
   + sub(i=-i,Qmm1m)*sft(Smp1m,-2) - sub(i=-i,Smm1m)*sft(Qmp1m,-2)
   + sub(i=-i,Qmm2m)*sft(Smm,  -2) - sub(i=-i,Smm2m)*sft(Qmm,  -2)
   + sub(i=-i,Qmm3m)*sft(Smm1m,-2) - sub(i=-i,Smm3m)*sft(Qmm1m,-2)
   + sub(i=-i,Qmm4m)*sft(Smm2m,-2) - sub(i=-i,Smm4m)*sft(Qmm2m,-2)$
Cmmm2 := coeffn(Cmmm2,eps,0) + eps*coeffn(Cmmm2,eps,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate m,m-3 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ammm3 :=
   + sub(i=-i,Pmm)  *sft(Rmp3m,-3) - sub(i=-i,Rmm)  *sft(Pmp3m,-3)
   + sub(i=-i,Pmp1m)*sft(Rmp4m,-3) - sub(i=-i,Rmp1m)*sft(Pmp4m,-3)
   + sub(i=-i,Pmm1m)*sft(Rmp2m,-3) - sub(i=-i,Rmm1m)*sft(Pmp2m,-3)
   + sub(i=-i,Pmm2m)*sft(Rmp1m,-3) - sub(i=-i,Rmm2m)*sft(Pmp1m,-3)
   + sub(i=-i,Pmm3m)*sft(Rmm,  -3) - sub(i=-i,Rmm3m)*sft(Pmm,  -3)
   + sub(i=-i,Pmm4m)*sft(Rmm1m,-3) - sub(i=-i,Rmm4m)*sft(Pmm1m,-3)$
Ammm3 := coeffn(Ammm3,eps,0) + eps*coeffn(Ammm3,eps,1);

Bmmm3 :=
   + sub(i=-i,Pmm)  *sft(Smp3m,-3) - sub(i=-i,Rmm)  *sft(Qmp3m,-3)
   + sub(i=-i,Pmp1m)*sft(Smp4m,-3) - sub(i=-i,Rmp1m)*sft(Qmp4m,-3)
   + sub(i=-i,Pmp2m)*sft(Smp2m,-3) - sub(i=-i,Rmp2m)*sft(Qmp2m,-3)
   + sub(i=-i,Pmm2m)*sft(Smp1m,-3) - sub(i=-i,Rmm2m)*sft(Qmp1m,-3)
   + sub(i=-i,Pmm3m)*sft(Smm,  -3) - sub(i=-i,Rmm3m)*sft(Qmm,  -3)
   + sub(i=-i,Pmm4m)*sft(Smm1m,-3) - sub(i=-i,Rmm4m)*sft(Qmm1m,-3)$
Bmmm3 := coeffn(Bmmm3,eps,0) + eps*coeffn(Bmmm3,eps,1);

Cmmm3 :=
   + sub(i=-i,Qmm)  *sft(Smp3m,-3) - sub(i=-i,Smm)  *sft(Qmp3m,-3)
   + sub(i=-i,Qmp1m)*sft(Smp4m,-3) - sub(i=-i,Smp1m)*sft(Qmp4m,-3)
   + sub(i=-i,Qmm1m)*sft(Smp2m,-3) - sub(i=-i,Smm1m)*sft(Qmp2m,-3)
   + sub(i=-i,Qmm2m)*sft(Smp1m,-3) - sub(i=-i,Smm2m)*sft(Qmp1m,-3)
   + sub(i=-i,Qmm3m)*sft(Smm,  -3) - sub(i=-i,Smm3m)*sft(Qmm,  -3)
   + sub(i=-i,Qmm4m)*sft(Smm1m,-3) - sub(i=-i,Smm4m)*sft(Qmm1m,-3)$
Cmmm3 := coeffn(Cmmm3,eps,0) + eps*coeffn(Cmmm3,eps,1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate special elements of solution vectors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P00 := P000_target + eps**2*P002_target$
Q00 := sub(m=0,Qmm0_target) + eps**2*sub(m=0,Qmm2_target)$
R00 := R000_target + eps**2*R002_target$
S00 := sub(m=0,Smm0_target) + eps**2*sub(m=0,Smm2_target);

P0p1 := eps * P00p11_target$
Q0p1 := eps * sub(m=0,Qmmp11_target)$
R0p1 := eps * R00p1_target$
S0p1 := eps * sub(m=0,Smmp1_target)$

P0p2 := eps * P00p21_target$
Q0p2 := eps * sub(m=0,Qmmp21_target)$
R0p2 := eps * R00p2_target$
S0p2 := eps * sub(m=0,Smmp2_target)$

P0p3 := eps * P00p31_target$
Q0p3 := eps * sub(m=0,Qmmp31_target)$
R0p3 := eps * R00p3_target$
S0p3 := eps * sub(m=0,Smmp3_target)$

P0p4 := eps * P00p41_target$
Q0p4 := eps * sub(m=0,Qmmp41_target)$
R0p4 := eps * R00p4_target$
S0p4 := eps * sub(m=0,Smmp4_target)$

P0m1 := eps * P00m11_target$
Q0m1 := eps * sub(m=0,Qmmm11_target)$
R0m1 := eps * R00m1_target$
S0m1 := eps * sub(m=0,Smmm1_target)$

P0m2 := eps * P00m21_target$
Q0m2 := eps * sub(m=0,Qmmm21_target)$
R0m2 := eps * R00m2_target$
S0m2 := eps * sub(m=0,Smmm2_target)$

P0m3 := eps * P00m31_target$
Q0m3 := eps * sub(m=0,Qmmm31_target)$
R0m3 := eps * R00m3_target$
S0m3 := eps * sub(m=0,Smmm3_target)$

P0m4 := eps * P00m41_target$
Q0m4 := eps * sub(m=0,Qmmm41_target)$
R0m4 := eps * R00m4_target$
S0m4 := eps * sub(m=0,Smmm4_target)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate 0,0 elements of coupling matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A00 :=
   + sub(i=-i,P00) *R00  - sub(i=-i,R00) *P00
   + sub(i=-i,P0p1)*R0p1 - sub(i=-i,R0p1)*P0p1
   + sub(i=-i,P0p2)*R0p2 - sub(i=-i,R0p2)*P0p2
   + sub(i=-i,P0p3)*R0p3 - sub(i=-i,R0p3)*P0p3
   + sub(i=-i,P0p4)*R0p4 - sub(i=-i,R0p4)*P0p4
   + sub(i=-i,P0m1)*R0m1 - sub(i=-i,R0m1)*P0m1
   + sub(i=-i,P0m2)*R0m2 - sub(i=-i,R0m2)*P0m2
   + sub(i=-i,P0m3)*R0m3 - sub(i=-i,R0m3)*P0m3
   + sub(i=-i,P0m4)*R0m4 - sub(i=-i,R0m4)*P0m4;

B00 :=
   + sub(i=-i,P00) *S00  - sub(i=-i,R00) *Q00
   + sub(i=-i,P0p1)*S0p1 - sub(i=-i,R0p1)*Q0p1
   + sub(i=-i,P0p2)*S0p2 - sub(i=-i,R0p2)*Q0p2
   + sub(i=-i,P0p3)*S0p3 - sub(i=-i,R0p3)*Q0p3
   + sub(i=-i,P0p4)*S0p4 - sub(i=-i,R0p4)*Q0p4
   + sub(i=-i,P0m1)*S0m1 - sub(i=-i,R0m1)*Q0m1
   + sub(i=-i,P0m2)*S0m2 - sub(i=-i,R0m2)*Q0m2
   + sub(i=-i,P0m3)*S0m3 - sub(i=-i,R0m3)*Q0m3
   + sub(i=-i,P0m4)*S0m4 - sub(i=-i,R0m4)*Q0m4;

C00 :=
   + sub(i=-i,Q00) *S00  - sub(i=-i,S00) *Q00
   + sub(i=-i,Q0p1)*S0p1 - sub(i=-i,S0p1)*Q0p1
   + sub(i=-i,Q0p2)*S0p2 - sub(i=-i,S0p2)*Q0p2
   + sub(i=-i,Q0p3)*S0p3 - sub(i=-i,S0p3)*Q0p3
   + sub(i=-i,Q0p4)*S0p4 - sub(i=-i,S0p4)*Q0p4
   + sub(i=-i,Q0m1)*S0m1 - sub(i=-i,S0m1)*Q0m1
   + sub(i=-i,Q0m2)*S0m2 - sub(i=-i,S0m2)*Q0m2
   + sub(i=-i,Q0m3)*S0m3 - sub(i=-i,S0m3)*Q0m3
   + sub(i=-i,Q0m4)*S0m4 - sub(i=-i,S0m4)*Q0m4;

;bye;
