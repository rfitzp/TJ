% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recreation of Jim Hastie and Richard Fitzpatrick's T7 calculation
% Extended for TJ calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out "TJ.out"$

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to evaluate parallel gradient
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure lgrd(ex)$ df(ex,t) - i*n*q*ex$

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fundamental equilibrium quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depend H1,r$
depend H2,r$
depend H3,r$
depend H4,r$
depend V2,r$
depend V3,r$
depend V4,r$
depend p2,r$
depend x,r,w$
depend z,r,w$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express x and z in terms of flux-surface label r, geometric
% angle w, and previously defined quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x := - r*cos(w) + eps*H1
       		+ eps*H2*cos(w) + eps*H3*cos(2*w) + eps*H4*cos(3*w)
		+ eps*V2*sin(w) + eps*V3*sin(2*w) + eps*V4*sin(3*w)
		+ eps**2*p2*cos(w)$
z :=   r*sin(w) + eps*H2*sin(w) + eps*H3*sin(2*w) + eps*H4*sin(3*w)
                - eps*V2*cos(w) - eps*V3*cos(2*w) - eps*V4*cos(3*w)
       		- eps**2*p2*sin(w)$

% %%%%%%%%%%%%%%%%%
% Evaluate Jacobian
% %%%%%%%%%%%%%%%%%
xr := df(x,r)$
xw := df(x,w)$
zr := df(z,r)$
zw := df(z,w)$
j  := xw*zr - xr*zw$

write "Jacobian in geometric coordinates:";
write "j := j0 + eps*j1 + eps**2*j2";
j0 := coeffn(j,eps,0);
j1 := coeffn(j,eps,1);
j2 := coeffn(j,eps,2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make transformation to straight field-line coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1  := 1 + eps*x$
ix1 := inv(x1,eps)$
k   := j*ix1$
ik  := int(k,w)$
ik0 := (sub(w=2*pi, ik) - sub(w=0,ik)) /(2*pi)$

write "j/x := k0 + eps*k1 + eps*eps*k2:"$
k0 := coeffn(k,eps,0);
k1 := coeffn(k,eps,1);
k2 := coeffn(k,eps,2);

write "Choice of P2 with associated residual:";
sp0  :=   (2-1)*H2**2/(2*r) + (3-1)*H3**2/(2*r) + (4-1)*H4**2/(2*r)
        + (2-1)*V2**2/(2*r) + (3-1)*V3**2/(2*r) + (4-1)*V4**2/(2*r)$
q2   := r**3/8 - r*H1/2 - sp0;
p2   := q2$
res2 := coeffn(ik0,eps,2);

write "Straight field-line flux-surface label:";
rst := sqrt(2) * sqrt(int(ik0,r));

write "theta := tt0 + eps*tt1 + eps**2*tt2";
theta := ik/r$
tt0   := coeffn(theta,eps,0);
tt1   := coeffn(theta,eps,1);
tt2   := coeffn(theta,eps,2);

write "Secular variation residual:";
rest := (sub(w=2*pi,theta) - sub(w=0,theta))/(2*pi) - 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate covariant metric tensor elements in geometic coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Evaluate covariant metric tensor in geometric coordinates:";
gww := zw**2 + xw**2$
grw := zr*zw + xr*xw$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate contravariant  metric tensor elements in geometic coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Evaluate contravariant metric tensor in geometric coordinates:";
ij     := inv(j,eps)$
ij2    := ij**2$
gradr2 := gww*ij2$
grrgrw := - grw*ij2$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate contravariant metric tensor elements in straight field-line coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Transform contravariant metric tensor to straight coordinates:";
tw     := df(theta,w)$
ttr    := df(theta,r)$
grrgrt := grrgrw*tw + gradr2*ttr$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate geometric angle w as function of straight field-line angle t
%
%  t  = w + eps*www1(w) + O(eps**2)
%  t  = w + eps*www1(t) + O(eps**2)
%  w  = t + eps*ttt1(t) + O(eps**2)
%
%  ttt1(t) = - www1(t)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Express geometric angle in terms of straight angle:";
www1 := coeffn(theta,eps,1)$
ttt1 := - sub(w=t,www1)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express all trigonometric functions in terms of t
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let           cos(w)   = cos(t)   -   sin(t)  *eps*ttt1$
for all n let cos(n*w) = cos(n*t) - n*sin(n*t)*eps*ttt1$
let           sin(w)   = sin(t)   +   cos(t)  *eps*ttt1$
for all n let sin(n*w) = sin(n*t) + n*cos(n*t)*eps*ttt1$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output metric tensor elements with t-averaged eps**2 cooefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Contravariant metric tensor elements in straight coordinates:";

write "gr2 := gr20 + eps*gr21 + eps**2*gr22";
gr2  := av2(gradr2,eps,t)$
gr20 := coeffn(gr2,eps,0);
gr21 := coeffn(gr2,eps,1);
gr22 := coeffn(gr2,eps,2);

gr20_target := 1$
gr21_target := 2*df(H1,r)*cos(t)
   + 2*df(H2,r)*cos(2*t) + 2*df(H3,r)*cos(3*t) + 2*df(H4,r)*cos(4*t)
   + 2*df(V2,r)*sin(2*t) + 2*df(V3,r)*sin(3*t) + 2*df(V4,r)*sin(4*t)$
gr22_target := 3*r*r/4 - H1
   + (df(H1,r)*df(H1,r) + (1*1-1)*H1*H1/r/r)/2
   + (df(H2,r)*df(H2,r) + (2*2-1)*H2*H2/r/r)/2
   + (df(H3,r)*df(H3,r) + (3*3-1)*H3*H3/r/r)/2
   + (df(H4,r)*df(H4,r) + (4*4-1)*H4*H4/r/r)/2
   + (df(V2,r)*df(V2,r) + (2*2-1)*V2*V2/r/r)/2
   + (df(V3,r)*df(V3,r) + (3*3-1)*V3*V3/r/r)/2
   + (df(V4,r)*df(V4,r) + (4*4-1)*V4*V4/r/r)/2$

write "Residuals of gr2 elements:";
gr20 - gr20_target;
gr21 - gr21_target;
gr22 - gr22_target;

write "grt2 := grt20 + eps*grt21 + eps**2*grt22";
grt2  := av2(grrgrt,eps,t)$
grt20 := coeffn(grt2,eps,0);
grt21 := coeffn(grt2,eps,1);
grt22 := coeffn(grt2,eps,2);
grt22 := 0;

grt20_target := 0$
grt21_target := sin(t)
   - (df(H1,r,2) + df(H1,r)/r + (1*1-1)*H1/r/r) * sin(t)/1
   - (df(H2,r,2) + df(H2,r)/r + (2*2-1)*H2/r/r) * sin(2*t)/2
   - (df(H3,r,2) + df(H3,r)/r + (3*3-1)*H3/r/r) * sin(3*t)/3
   - (df(H4,r,2) + df(H4,r)/r + (4*4-1)*H4/r/r) * sin(4*t)/4
   + (df(V2,r,2) + df(V2,r)/r + (2*2-1)*V2/r/r) * cos(2*t)/2
   + (df(V3,r,2) + df(V3,r)/r + (3*3-1)*V3/r/r) * cos(3*t)/3
   + (df(V4,r,2) + df(V4,r)/r + (4*4-1)*V4/r/r) * cos(4*t)/4$
grt22_target := 0$

write "Residuals of grt2 elements:";
grt20 - grt20_target;
grt21 - grt21_target;
grt22 - grt22_target;

write "x*x := x20 + eps*x21 + eps**2*x22";
x1  := 1 + eps*x$
x2  := x1*x1$
x2  := av2(x2,eps,t)$
x20 := coeffn(x2,eps,0);
x21 := coeffn(x2,eps,1);
x22 := coeffn(x2,eps,2);

x20_target := 1$
x21_target := - 2 * r * cos(t)$
x22_target := - r*r/2 + r*df(H1,r) + 2*H1$

write "Residuals of x2 elements:";
x20 - x20_target;
x21 - x21_target;
x22 - x22_target;

write "x2gr2 := x2gr20 + eps*x2gr21 + eps**2*x2gr22";
x2gr2  := av2(x2*gr2,eps,t)$
x2gr20 := coeffn(x2gr2,eps,0);
x2gr21 := coeffn(x2gr2,eps,1);
x2gr22 := coeffn(x2gr2,eps,2);

x2gr20_target := 1$
x2gr21_target := - 2*r*cos(t)
   + 2*df(H1,r) * cos(t)
   + 2*df(H2,r) * cos(2*t)
   + 2*df(H3,r) * cos(3*t)
   + 2*df(H4,r) * cos(4*t)
   + 2*df(V2,r) * sin(2*t)
   + 2*df(V3,r) * sin(3*t)
   + 2*df(V4,r) * sin(4*t)$
x2gr22_target := r*r/4 - df(H1,r)*r + H1
   + (df(H1,r)*df(H1,r) + (1*1-1)*H1*H1/r/r)/2
   + (df(H2,r)*df(H2,r) + (2*2-1)*H2*H2/r/r)/2
   + (df(H3,r)*df(H3,r) + (3*3-1)*H3*H3/r/r)/2
   + (df(H4,r)*df(H4,r) + (4*4-1)*H4*H4/r/r)/2
   + (df(V2,r)*df(V2,r) + (2*2-1)*V2*V2/r/r)/2
   + (df(V3,r)*df(V3,r) + (3*3-1)*V3*V3/r/r)/2
   + (df(V4,r)*df(V4,r) + (4*4-1)*V4*V4/r/r)/2$

write "Residuals of x2gr2 elements:";
x2gr20 - x2gr20_target;
x2gr21 - x2gr21_target;
x2gr22 - x2gr22_target;

write "x2grt2  := x2grt20 + eps*x2grt21 + eps**2*x2grt22";
x2grt2  := av2(x2*grt2,eps,t)$
x2grt20 := coeffn(x2grt2,eps,0);
x2grt21 := coeffn(x2grt2,eps,1);
x2grt22 := coeffn(x2grt2,eps,2);
x2grt22 := 0;

x2grt20_target := 0$
x2grt21_target := sin(t)
   - (df(H1,r,2) + df(H1,r)/r + (1*1-1)*H1/r/r) * sin(t)/1
   - (df(H2,r,2) + df(H2,r)/r + (2*2-1)*H2/r/r) * sin(2*t)/2
   - (df(H3,r,2) + df(H3,r)/r + (3*3-1)*H3/r/r) * sin(3*t)/3
   - (df(H4,r,2) + df(H4,r)/r + (4*4-1)*H4/r/r) * sin(4*t)/4
   + (df(V2,r,2) + df(V2,r)/r + (2*2-1)*V2/r/r) * cos(2*t)/2
   + (df(V3,r,2) + df(V3,r)/r + (3*3-1)*V3/r/r) * cos(3*t)/3
   + (df(V4,r,2) + df(V4,r)/r + (4*4-1)*V4/r/r) * cos(4*t)/4$
x2grt22_target := 0$

write "Residuals of x2grt2 elements:";
x2grt20 - x2grt20_target;
x2grt21 - x2grt21_target;
x2grt22 - x2grt22_target;

igr2 := av2(inv(gr2,eps),eps,t)$
ix2  := av2(inv(x2,eps),eps,t)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expand Grad-Shafranov equation to order eps**2
% Plasma pressure and safety-factor are not expanded
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depend q,r$
depend p,r$
depend g2,r$
depend g4,r$

gg  := 1 + eps**2*g2$
ggp := df(g2,r) + eps**2*df(g4,r)$
pp  := df(p,r)$

write "Express Grad-Shafranov equation in straight coordinates:";
fac1 := r**2*gg*gr2/q$
fac2 := r**2*gg*grt2/q$
gs   := (gg/q) * (df(fac1,r) + df(fac2,t)) + x2*pp + gg*ggp$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier analyze Grad-Shafranov equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Fourier analyze Grad-Shafranov equation:";
gs0  := fcos(gs,t,0)$
gs00 := coeffn(gs0,eps,0)$
gs02 := coeffn(gs0,eps,2)$
gs1c := coeffn(fcos(gs,t,1),eps,1)$
gs2c := coeffn(fcos(gs,t,2),eps,1)$
gs3c := coeffn(fcos(gs,t,3),eps,1)$
gs4c := coeffn(fcos(gs,t,4),eps,1)$
gs1s := coeffn(fsin(gs,t,1),eps,1)$
gs2s := coeffn(fsin(gs,t,2),eps,1)$
gs3s := coeffn(fsin(gs,t,3),eps,1)$
gs4s := coeffn(fsin(gs,t,4),eps,1)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output components of Grad-Shafranov equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Components of Grad-Shafranov equation:";

write "Zeroth-order equilibrium equation:";
c0   := coeffn(gs00,df(g2,r),0);
c1   := coeffn(gs00,df(g2,r),1);
g2p  := -c0/c1;

write "l=1 shaping equation:";
c0   := coeffn(gs1c,df(H1,r,2),0)$
c1   := coeffn(gs1c,df(H1,r,2),1)$
H1pp := -c0/c1;

write "l=2 shaping equation:";
c0   := coeffn(gs2c,df(H2,r,2),0)$
c1   := coeffn(gs2c,df(H2,r,2),1)$
H2pp := -c0/c1;
c0   := coeffn(gs2s,df(V2,r,2),0)$
c1   := coeffn(gs2s,df(V2,r,2),1)$
V2pp := -c0/c1;

write "l=3 shaping equation:";
c0   := coeffn(gs3c,df(H3,r,2),0)$
c1   := coeffn(gs3c,df(H3,r,2),1)$
H3pp := -c0/c1;
c0   := coeffn(gs3s,df(V3,r,2),0)$
c1   := coeffn(gs3s,df(V3,r,2),1)$
V3pp := -c0/c1;

write "l=4 shaping equation:";
c0   := coeffn(gs4c,df(H4,r,2),0)$
c1   := coeffn(gs4c,df(H4,r,2),1)$
H4pp := -c0/c1;
c0   := coeffn(gs4s,df(V4,r,2),0)$
c1   := coeffn(gs4s,df(V4,r,2),1)$
V4pp := -c0/c1;

write "Second-order equilibrium equation:";
c0   := coeffn(gs02,df(g4,r),0)$
c1   := coeffn(gs02,df(g4,r),1)$
g4p  := -c0/c1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only need to retain harmonic terms up to l=4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let cos(5*t) = 0$
let cos(6*t) = 0$
let cos(7*t) = 0$
let cos(8*t) = 0$
let sin(5*t) = 0$
let sin(6*t) = 0$
let sin(7*t) = 0$
let sin(8*t) = 0$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to eliminate harmonic terms in second order component of
% polynomial ex in eps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure qav2(ex)$
   begin scalar c0,c1,c2$
      c0 := coeffn(ex,eps,0)$
      c1 := coeffn(ex,eps,1)$
      c2 := coeffn(ex,eps,2)$
      c2 := sub(cos(t)=0,cos(2*t)=0,cos(3*t)=0,cos(4*t)=0,
	        sin(t)=0,sin(2*t)=0,sin(3*t)=0,sin(4*t)=0,
	        c2)$
      return c0 + eps*c1 + eps*eps*c2$
end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to eliminate harmonic terms in expression ex
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure av(ex)$
   begin scalar c0$
      c0 := ex$
      c0 := sub(cos(t)=0,cos(2*t)=0,cos(3*t)=0,cos(4*t)=0,
	       sin(t)=0,sin(2*t)=0,sin(3*t)=0,sin(4*t)=0,
	       c0)$
      return c0$
end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate factors appearing in outer region pdes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depend g,r$
depend invf,r$
gp := df(g,r)$
gf := gp*invf$
pf := pp*(invf**2)$

write "S/r := ss";
ss := i*n*r;

write "Q/r := qq0 + eps*qq1 + eps**2*qq2";
qq  := igr2/(i*n*r)$
qq0 := coeffn(qq,eps,0);
qq1 := coeffn(qq,eps,1);
qq2 := coeffn(qq,eps,2);

write "T/r := tc0 + eps*tc1 + eps**2*tc2";
tt  := grt2*igr2 - qq*gf$
tt  := qav2(tt)$
tc0 := coeffn(tt,eps,0);
tc1 := coeffn(tt,eps,1);
tc2 := coeffn(tt,eps,2);
ts  := sub(i=-i,tt)$

write "U/r := uu0 + eps*uu1 + eps**2*uu2";
uu  := igr2*x2*pf$
uu  := qav2(uu)$
uu0 := coeffn(uu,eps,0);
uu1 := coeffn(uu,eps,1);
uu2 := coeffn(uu,eps,2);

write "V/r := vv0 + eps*vv1 + eps**2*vv2";
vv  := qq * (gf**2 - n**2*ix2)$
vv  := qav2(vv)$
vv0 := coeffn(vv,eps,0);
vv1 := coeffn(vv,eps,1);
vv2 := coeffn(vv,eps,2);

write "W/r := ww0 + eps*ww1 + eps**2*ww2";
ww  := 2*gf*uu - df(gf,r)$
ww  := qav2(ww)$
ww0 := coeffn(ww,eps,0);
ww1 := coeffn(ww,eps,1);
ww2 := coeffn(ww,eps,2);

write "X/r := xx0 + eps*xx1 + eps**2*xx2";
tx  := ts*x2$
fr  := r*invf$
xx  := ss*pf * (df(tx,t) + df(x2,r) + x2*df(fr,r)/fr - x2*uu)$
xx  := qav2(xx)$
xx0 := coeffn(xx,eps,0);
xx1 := coeffn(xx,eps,1);
xx2 := coeffn(xx,eps,2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate primitive coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arg := exp( i*(m+jj)*t)$
mul := exp(-i*m     *t)$

bmj := qq*df(arg,t)$
bmj := df(bmj,t) + eps**2*ss*arg$
bmj := mul*bmj$

cmj := tt*lgrd(arg) + uu*arg$
cmj := - df(cmj,t)$
cmj := mul*cmj$

dmj := ts*df(arg,t)$
dmj := - lgrd(dmj) + uu*df(arg,t)$
dmj := mul*dmj$

emj := vv*lgrd(arg)$
emj := - lgrd(emj) + ww*lgrd(arg) + xx*arg$
emj := mul*emj$

for all x,y let e**(x+y)   = e**x * e**y$
let             e**(i*t)   = cos(t)   + i*sin(t)$
for all n let   e**(i*n*t) = cos(n*t) + i*sin(n*t)$

off revpri$
factor m$;
order m,n,g,df(g,r),df(g,r,2),df(p,r),df(p,r,2),invf,df(invf,r),df(invf,r,2),q,r$
order H1,df(H1,r),df(H1,r,2),
      H2,df(H2,r),df(H2,r,2),V2,df(V2,r),df(V2,r,2),
      H3,df(H3,r),df(H3,r,2),V3,df(V3,r),df(V3,r,2),
      H4,df(H4,r),df(H4,r,2),V4,df(V4,r),df(V4,r,2)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=0 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = 0$
let q = (m - nmq0)/n$
factor nmq0$
bm0 := -i*r*bmj$
cm0 := -i*r*cmj$
dm0 := -i*r*dmj$
em0 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = 0:";

write "bm0 := bm00 + eps*bm01 + eps**2*bm02";
bm00 := coeffn(bm0,eps,0);
bm01 := coeffn(bm0,eps,1);
bm02 := coeffn(bm0,eps,2);

write "cm0 := cm00 + eps*cm01 + eps**2*cm02";
cm00 := coeffn(cm0,eps,0);
cm01 := coeffn(cm0,eps,1);
cm02 := coeffn(cm0,eps,2);

write "dm0 := dm00 + eps*dm01 + eps**2*dm02";
dm00 := coeffn(dm0,eps,0);
dm01 := coeffn(dm0,eps,1);
dm02 := coeffn(dm0,eps,2);

write "em0 := em00 + eps*em01 + eps**2*em02";
em00 := coeffn(em0,eps,0);
em01 := coeffn(em0,eps,1);
em02 := coeffn(em0,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = 1$
clear q$
let q = (m + 1 - nmqp1)/n$
factor nmqp1$
bmp1 := -i*r*bmj$
cmp1 := -i*r*cmj$
dmp1 := -i*r*dmj$
emp1 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = +1:";

write "bmp1 := bmp10 + eps*bp101 + eps**2*bmp12";
bmp10 := coeffn(bmp1,eps,0);
bmp11 := coeffn(bmp1,eps,1);
bmp12 := coeffn(bmp1,eps,2);

write "cmp1 := cmp10 + eps*cmp11 + eps**2*cmp12";
cmp10 := coeffn(cmp1,eps,0);
cmp11 := coeffn(cmp1,eps,1);
cmp12 := coeffn(cmp1,eps,2);

write "dmp1 := dmp10 + eps*dmp11 + eps**2*dmp12";
dmp10 := coeffn(dmp1,eps,0);
dmp11 := coeffn(dmp1,eps,1);
dmp12 := coeffn(dmp1,eps,2);

write "emp1 := emp10 + eps*emp11 + eps**2*emp12";
emp10 := coeffn(emp1,eps,0);
emp11 := coeffn(emp1,eps,1);
emp12 := coeffn(emp1,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = -1$
clear q$
let q = (m - 1 - nmqm1)/n$
factor nmqm1$
bmm1 := -i*r*bmj$
cmm1 := -i*r*cmj$
dmm1 := -i*r*dmj$
emm1 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = -1:";

write "bmm1 := bmm10 + eps*bm101 + eps**2*bmm12";
bmm10 := coeffn(bmm1,eps,0);
bmm11 := coeffn(bmm1,eps,1);
bmm12 := coeffn(bmm1,eps,2);

write "cmm1 := cmm10 + eps*cmm11 + eps**2*cmm12";
cmm10 := coeffn(cmm1,eps,0);
cmm11 := coeffn(cmm1,eps,1);
cmm12 := coeffn(cmm1,eps,2);

write "dmm1 := dmm10 + eps*dmm11 + eps**2*dmm12";
dmm10 := coeffn(dmm1,eps,0);
dmm11 := coeffn(dmm1,eps,1);
dmm12 := coeffn(dmm1,eps,2);

write "emm1 := emm10 + eps*emm11 + eps**2*emm12";
emm10 := coeffn(emm1,eps,0);
emm11 := coeffn(emm1,eps,1);
emm12 := coeffn(emm1,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = 2$
clear q$
let q = (m + 2 - nmqp2)/n$
factor nmqp2$
bmp2 := -i*r*bmj$
cmp2 := -i*r*cmj$
dmp2 := -i*r*dmj$
emp2 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = +2:";

write "bmp2 := bmp20 + eps*bp201 + eps**2*bmp22";
bmp20 := coeffn(bmp2,eps,0);
bmp21 := coeffn(bmp2,eps,1);
bmp22 := coeffn(bmp2,eps,2);

write "cmp2 := cmp20 + eps*cmp21 + eps**2*cmp22";
cmp20 := coeffn(cmp2,eps,0);
cmp21 := coeffn(cmp2,eps,1);
cmp22 := coeffn(cmp2,eps,2);

write "dmp2 := dmp20 + eps*dmp21 + eps**2*dmp22";
dmp20 := coeffn(dmp2,eps,0);
dmp21 := coeffn(dmp2,eps,1);
dmp22 := coeffn(dmp2,eps,2);

write "emp2 := emp20 + eps*emp21 + eps**2*emp22";
emp20 := coeffn(emp2,eps,0);
emp21 := coeffn(emp2,eps,1);
emp22 := coeffn(emp2,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = -2$
clear q$
let q = (m - 2 - nmqm2)/n$
factor nmqm2$
bmm2 := -i*r*bmj$
cmm2 := -i*r*cmj$
dmm2 := -i*r*dmj$
emm2 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = -2:";

write "bmm2 := bmm20 + eps*bm201 + eps**2*bmm22";
bmm20 := coeffn(bmm2,eps,0);
bmm21 := coeffn(bmm2,eps,1);
bmm22 := coeffn(bmm2,eps,2);

write "cmm2 := cmm20 + eps*cmm21 + eps**2*cmm22";
cmm20 := coeffn(cmm2,eps,0);
cmm21 := coeffn(cmm2,eps,1);
cmm22 := coeffn(cmm2,eps,2);

write "dmm2 := dmm20 + eps*dmm21 + eps**2*dmm22";
dmm20 := coeffn(dmm2,eps,0);
dmm21 := coeffn(dmm2,eps,1);
dmm22 := coeffn(dmm2,eps,2);

write "emm2 := emm20 + eps*emm21 + eps**2*emm22";
emm20 := coeffn(emm2,eps,0);
emm21 := coeffn(emm2,eps,1);
emm22 := coeffn(emm2,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = 3$
clear q$
let q = (m + 3 - nmqp3)/n$
factor nmqp3$
bmp3 := -i*r*bmj$
cmp3 := -i*r*cmj$
dmp3 := -i*r*dmj$
emp3 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = +3:";

write "bmp3 := bmp30 + eps*bp301 + eps**2*bmp32";
bmp30 := coeffn(bmp3,eps,0);
bmp31 := coeffn(bmp3,eps,1);
bmp32 := coeffn(bmp3,eps,2);

write "cmp3 := cmp30 + eps*cmp31 + eps**2*cmp32";
cmp30 := coeffn(cmp3,eps,0);
cmp31 := coeffn(cmp3,eps,1);
cmp32 := coeffn(cmp3,eps,2);

write "dmp3 := dmp30 + eps*dmp31 + eps**2*dmp32";
dmp30 := coeffn(dmp3,eps,0);
dmp31 := coeffn(dmp3,eps,1);
dmp32 := coeffn(dmp3,eps,2);

write "emp3 := emp30 + eps*emp31 + eps**2*emp32";
emp30 := coeffn(emp3,eps,0);
emp31 := coeffn(emp3,eps,1);
emp32 := coeffn(emp3,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = -3$
clear q$
let q = (m - 3 - nmqm3)/n$
factor nmqm3$
bmm3 := -i*r*bmj$
cmm3 := -i*r*cmj$
dmm3 := -i*r*dmj$
emm3 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = -3:";

write "bmm3 := bmm30 + eps*bm301 + eps**2*bmm32";
bmm30 := coeffn(bmm3,eps,0);
bmm31 := coeffn(bmm3,eps,1);
bmm32 := coeffn(bmm3,eps,2);

write "cmm3 := cmm30 + eps*cmm31 + eps**2*cmm32";
cmm30 := coeffn(cmm3,eps,0);
cmm31 := coeffn(cmm3,eps,1);
cmm32 := coeffn(cmm3,eps,2);

write "dmm3 := dmm30 + eps*dmm31 + eps**2*dmm32";
dmm30 := coeffn(dmm3,eps,0);
dmm31 := coeffn(dmm3,eps,1);
dmm32 := coeffn(dmm3,eps,2);

write "emm3 := emm30 + eps*emm31 + eps**2*emm32";
emm30 := coeffn(emm3,eps,0);
emm31 := coeffn(emm3,eps,1);
emm32 := coeffn(emm3,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = 4$
clear q$
let q = (m + 4 - nmqp4)/n$
factor nmqp4$
bmp4 := -i*r*bmj$
cmp4 := -i*r*cmj$
dmp4 := -i*r*dmj$
emp4 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = +4:";

write "bmp4 := bmp40 + eps*bp401 + eps**2*bmp42";
bmp40 := coeffn(bmp4,eps,0);
bmp41 := coeffn(bmp4,eps,1);
bmp42 := coeffn(bmp4,eps,2);

write "cmp4 := cmp40 + eps*cmp41 + eps**2*cmp42";
cmp40 := coeffn(cmp4,eps,0);
cmp41 := coeffn(cmp4,eps,1);
cmp42 := coeffn(cmp4,eps,2);

write "dmp4 := dmp40 + eps*dmp41 + eps**2*dmp42";
dmp40 := coeffn(dmp4,eps,0);
dmp41 := coeffn(dmp4,eps,1);
dmp42 := coeffn(dmp4,eps,2);

write "emp4 := emp40 + eps*emp41 + eps**2*emp42";
emp40 := coeffn(emp4,eps,0);
emp41 := coeffn(emp4,eps,1);
emp42 := coeffn(emp4,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
let jj = -4$
clear q$
let q = (m - 4 - nmqm4)/n$
factor nmqm4$
bmm4 := -i*r*bmj$
cmm4 := -i*r*cmj$
dmm4 := -i*r*dmj$
emm4 := -i*r*emj$
for all x let cos(x) = 0$
for all x let sin(x) = 0$

write "Primitive coupling coefficients: j = -4:";

write "bmm4 := bmm40 + eps*bm401 + eps**2*bmm42";
bmm40 := coeffn(bmm4,eps,0);
bmm41 := coeffn(bmm4,eps,1);
bmm42 := coeffn(bmm4,eps,2);

write "cmm4 := cmm40 + eps*cmm41 + eps**2*cmm42";
cmm40 := coeffn(cmm4,eps,0);
cmm41 := coeffn(cmm4,eps,1);
cmm42 := coeffn(cmm4,eps,2);

write "dmm4 := dmm40 + eps*dmm41 + eps**2*dmm42";
dmm40 := coeffn(dmm4,eps,0);
dmm41 := coeffn(dmm4,eps,1);
dmm42 := coeffn(dmm4,eps,2);

write "emm4 := emm40 + eps*emm41 + eps**2*emm42";
emm40 := coeffn(emm4,eps,0);
emm41 := coeffn(emm4,eps,1);
emm42 := coeffn(emm4,eps,2);

for all x clear cos(x)$
for all x clear sin(x)$

clear q$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform to final coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invf       := q*inv(gg,eps)/r$
df(g,r)    := ggp$
df(g2,r)   := g2p$
df(g4,r)   := g4p$
df(H1,r,2) := H1pp$
df(H2,r,2) := H2pp$
df(H3,r,2) := H3pp$
df(H4,r,2) := H4pp$
df(V2,r,2) := V2pp$
df(V3,r,2) := V3pp$
df(V4,r,2) := V4pp$

nmq0  := m     - n*q$
nmqp1 := m + 1 - n*q$
nmqm1 := m - 1 - n*q$
nmqp2 := m + 2 - n*q$
nmqm2 := m - 2 - n*q$
nmqp3 := m + 3 - n*q$
nmqm3 := m - 3 - n*q$
nmqp4 := m + 4 - n*q$
nmqm4 := m - 4 - n*q$

lam00 := - cm00/bm00$
lam02 := - (coeffn(cm02,i,0) - cm00*bm02/bm00)/bm00$
lam0  := lam00 + eps*eps*lam02$

lamp1 := sub(m=kk+1,lam00)$
lamp1 := sub(kk=m,  lamp1)$
lamm1 := sub(m=kk-1,lam00)$
lamm1 := sub(kk=m,  lamm1)$
lamp2 := sub(m=kk+2,lam00)$
lamp2 := sub(kk=m,  lamp2)$
lamm2 := sub(m=kk-2,lam00)$
lamm2 := sub(kk=m,  lamm2)$
lamp3 := sub(m=kk+3,lam00)$
lamp3 := sub(kk=m,  lamp3)$
lamm3 := sub(m=kk-3,lam00)$
lamm3 := sub(kk=m,  lamm3)$
lamp4 := sub(m=kk+4,lam00)$
lamp4 := sub(kk=m,  lamp4)$
lamm4 := sub(m=kk-4,lam00)$
lamm4 := sub(kk=m,  lamm4)$

% %%%%%%%%%%%%%
% Calculate kmp
% %%%%%%%%%%%%%
kmp  := lam0/n + gf/n$
kmp0 := coeffn(kmp,eps,0)$
kmp1 := coeffn(kmp,eps,1)$
kmp2 := coeffn(kmp,eps,2)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=0 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lm0  := av(n*bm0)$
lm00 := coeffn(lm0,eps,0)$
lm02 := coeffn(lm0,eps,2)$

sp1  :=   (3/2)*df(H1,r)**2
        + (3/2)*df(H2,r)**2 - ((2*2-1)/2)*(H2/r)**2
	+ (3/2)*df(H3,r)**2 - ((3*3-1)/2)*(H3/r)**2
	+ (3/2)*df(H4,r)**2 - ((4*4-1)/2)*(H4/r)**2
	+ (3/2)*df(V2,r)**2 - ((2*2-1)/2)*(V2/r)**2
	+ (3/2)*df(V3,r)**2 - ((3*3-1)/2)*(V3/r)**2
	+ (3/2)*df(V4,r)**2 - ((4*4-1)/2)*(V4/r)**2$

sp2  :=   2*(2*2-1) * (df(H2,r)**2 - (11/3)*df(H2,r)*H2/r + 2*2*(H2/r)**2 - (1-sh)*(df(H2,r)*H2/r + (1/3)*(H2/r)**2))
        + 2*(3*3-1) * (df(H3,r)**2 - (11/3)*df(H3,r)*H3/r + 3*3*(H3/r)**2 - (1-sh)*(df(H3,r)*H3/r + (1/3)*(H3/r)**2))
	+ 2*(4*4-1) * (df(H4,r)**2 - (11/3)*df(H4,r)*H4/r + 4*4*(H4/r)**2 - (1-sh)*(df(H4,r)*H4/r + (1/3)*(H4/r)**2))
        + 2*(2*2-1) * (df(V2,r)**2 - (11/3)*df(V2,r)*V2/r + 2*2*(V2/r)**2 - (1-sh)*(df(V2,r)*V2/r + (1/3)*(V2/r)**2))
        + 2*(3*3-1) * (df(V3,r)**2 - (11/3)*df(V3,r)*V3/r + 3*3*(V3/r)**2 - (1-sh)*(df(V3,r)*V3/r + (1/3)*(V3/r)**2))
	+ 2*(4*4-1) * (df(V4,r)**2 - (11/3)*df(V4,r)*V4/r + 4*4*(V4/r)**2 - (1-sh)*(df(V4,r)*V4/r + (1/3)*(V4/r)**2))$

sp4  :=   (2-1) * (2+sh/2) * (df(H2,r)*V2 - df(V2,r)*H2)/r
        + (3-1) * (2+sh/3) * (df(H3,r)*V3 - df(V3,r)*H3)/r
	+ (4-1) * (2+sh/4) * (df(H4,r)*V4 - df(V4,r)*H4)/r$

lm00_target := m*m$
lm02_target := m*m*(- 3*r*r/4 + H1 + sp1) + n*n*r*r$

mm0  := av(cm0 + lam0*bm0)$
mm00 := coeffn(mm0,eps,0)$
mm02 := coeffn(mm0,eps,2)$

mm00_target := 0$
mm02_target := i*m*nmq0*sp4$ 

nm0  := av(dm0 - lam0*bm0)$
nm00 := coeffn(nm0,eps,0)$
nm02 := coeffn(nm0,eps,2)$

nm00_target := 0$
nm02_target := i*m*nmq0*sp4$

pm0  := av(em0 + lam0*(dm0 - mm0 - r*n*df(q,r)) - (m-n*q)*r*df(lam0,r))/n$
pm00 := coeffn(pm0,eps,0)$
pm02 := coeffn(pm0,eps,2)$

sh   := r*df(q,r)/q$
r2q  := r**2/q$
jcy  := df(r2q,r)/r$

pm00_target := nmq0**2 + (nmq0/m)*(q*r)*df(jcy,r)$

res1 := (r**2/q**3)*(2-sh) + (sh/q)*(3*r**2/4 - H1 - sp1) - (2/q)*((3/2)*r**2 - H1 - df(H1,r)*r - (2/3)*sp1)$
res2 := (n/m)*r*df(r*df(r2q,r),r) - df(r2q,r)**2 - r*df(r*pp,r) + m**2*(7*r**2/4 - H1 - 3*r*df(H1,r) + sp1)$

pm02_target :=   (2*pp*r)*(1-q**2)
	       - (nmq0/m)*((2*pp*r)*(2-sh) + (r*q)*df(res1,r) - sp2)
	       + (nmq0/m)**2*res2$

write "Residuals of j=0 coupling coefficients:";
dlm00 := lm00 - lm00_target;
dlm02 := lm02 - lm02_target;
dmm00 := mm00 - mm00_target;
dmm02 := mm02 - mm02_target;
dnm00 := nm00 - nm00_target;
dnm02 := nm02 - nm02_target;
dpm00 := pm00 - pm00_target;
dpm02 := pm02 - pm02_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmp11 := av(n*bmp11)$
mmp11 := av(cmp11 + lamp1*bmp11)$
nmp11 := av(dmp11 - lam00*bmp11)$
pmp11 := av(emp11 + lamp1*dmp11 - lam00*mmp11)/n$

lmp11 := coeffn(lmp11,eps,0)$
mmp11 := coeffn(mmp11,eps,0)$
nmp11 := coeffn(nmp11,eps,0)$
pmp11 := coeffn(pmp11,eps,0)$

lmp11_target := - m*(m+1)*df(H1,r)$
mmp11_target := - m*nmq0*pp*q**2      + m*nmqp1   *(r + df(H1,r)*(1-sh))$
nmp11_target := - (m+1)*nmqp1*pp*q**2 + (m+1)*nmq0*(r + df(H1,r)*(1-sh))$
pmp11_target := - (1+sh)*pp*q**2      + nmq0*nmqp1*(r - df(H1,r))$

write "Residuals of j=+1 coupling coefficients:";
dlmp11 := lmp11 - lmp11_target;
dmmp11 := mmp11 - mmp11_target;
dnmp11 := nmp11 - nmp11_target;
dpmp11 := pmp11 - pmp11_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmm11 := av(n*bmm11)$
mmm11 := av(cmm11 + lamm1*bmm11)$
nmm11 := av(dmm11 - lam00*bmm11)$
pmm11 := av(emm11 + lamm1*dmm11 - lam00*mmm11)/n$

lmm11 := coeffn(lmm11,eps,0)$
mmm11 := coeffn(mmm11,eps,0)$
nmm11 := coeffn(nmm11,eps,0)$
pmm11 := coeffn(pmm11,eps,0)$

lmm11_target := - m*(m-1)*df(H1,r)$
mmm11_target :=   m*nmq0*pp*q**2      - m*nmqm1   *(r + df(H1,r)*(1-sh))$
nmm11_target :=   (m-1)*nmqm1*pp*q**2 - (m-1)*nmq0*(r + df(H1,r)*(1-sh))$
pmm11_target := - (1+sh)*pp*q**2      + nmq0*nmqm1*(r - df(H1,r))$

write "Residuals of j=-1 coupling coefficients:";
dlmm11 := lmm11 - lmm11_target;
dmmm11 := mmm11 - mmm11_target;
dnmm11 := nmm11 - nmm11_target;
dpmm11 := pmm11 - pmm11_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmp21 := av(n*bmp21)$
mmp21 := av(cmp21 + lamp2*bmp21)$
nmp21 := av(dmp21 - lam00*bmp21)$
pmp21 := av(emp21 + lamp2*dmp21 - lam00*mmp21)/n$

lmp21 := coeffn(lmp21,eps,0)$
mmp21 := coeffn(mmp21,eps,0)$
nmp21 := coeffn(nmp21,eps,0)$
pmp21 := coeffn(pmp21,eps,0)$

lmp21_target := - m*(m+2)*(df(H2,r)+i*df(V2,r))$
mmp21_target :=   (1/2)*m    *nmqp2*((1-sh)*(df(H2,r)+i*df(V2,r)) - (2*2-1)*(H2+i*V2)/r)$
nmp21_target :=   (1/2)*(m+2)*nmq0 *((1-sh)*(df(H2,r)+i*df(V2,r)) - (2*2-1)*(H2+i*V2)/r)$
pmp21_target := - nmqp2*nmq0*(df(H2,r)+i*df(V2,r))$

write "Residuals of j=+2 coupling coefficients:";
dlmp21 := lmp21 - lmp21_target;
dmmp21 := mmp21 - mmp21_target;
dnmp21 := nmp21 - nmp21_target;
dpmp21 := pmp21 - pmp21_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmm21 := av(n*bmm21)$
mmm21 := av(cmm21 + lamm2*bmm21)$
nmm21 := av(dmm21 - lam00*bmm21)$
pmm21 := av(emm21 + lamm2*dmm21 - lam00*mmm21)/n$

lmm21 := coeffn(lmm21,eps,0)$
mmm21 := coeffn(mmm21,eps,0)$
nmm21 := coeffn(nmm21,eps,0)$
pmm21 := coeffn(pmm21,eps,0)$

lmm21_target := - m*(m-2)*(df(H2,r)-i*df(V2,r))$
mmm21_target := - (1/2)*m    *nmqm2*((1-sh)*(df(H2,r)-i*df(V2,r)) - (2*2-1)*(H2-i*V2)/r)$
nmm21_target := - (1/2)*(m-2)*nmq0 *((1-sh)*(df(H2,r)-i*df(V2,r)) - (2*2-1)*(H2-i*V2)/r)$
pmm21_target := - nmqm2*nmq0*(df(H2,r)-i*df(V2,r))$

write "Residuals of j=-2 coupling coefficients:";
dlmm21 := lmm21 - lmm21_target;
dmmm21 := mmm21 - mmm21_target;
dnmm21 := nmm21 - nmm21_target;
dpmm21 := pmm21 - pmm21_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmp31 := av(n*bmp31)$
mmp31 := av(cmp31 + lamp3*bmp31)$
nmp31 := av(dmp31 - lam00*bmp31)$
pmp31 := av(emp31 + lamp3*dmp31 - lam00*mmp31)/n$

lmp31 := coeffn(lmp31,eps,0)$
mmp31 := coeffn(mmp31,eps,0)$
nmp31 := coeffn(nmp31,eps,0)$
pmp31 := coeffn(pmp31,eps,0)$

lmp31_target := - m*(m+3)*(df(H3,r)+i*df(V3,r))$
mmp31_target :=   (1/3)*m    *nmqp3*((1-sh)*(df(H3,r)+i*df(V3,r)) - (3*3-1)*(H3+i*V3)/r)$
nmp31_target :=   (1/3)*(m+3)*nmq0 *((1-sh)*(df(H3,r)+i*df(V3,r)) - (3*3-1)*(H3+i*V3)/r)$
pmp31_target := - nmqp3*nmq0*(df(H3,r)+i*df(V3,r))$

write "Residuals of j=+3 coupling coefficients:";
dlmp31 := lmp31 - lmp31_target;
dmmp31 := mmp31 - mmp31_target;
dnmp31 := nmp31 - nmp31_target;
dpmp31 := pmp31 - pmp31_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmm31 := av(n*bmm31)$
mmm31 := av(cmm31 + lamm3*bmm31)$
nmm31 := av(dmm31 - lam00*bmm31)$
pmm31 := av(emm31 + lamm3*dmm31 - lam00*mmm31)/n$

lmm31 := coeffn(lmm31,eps,0)$
mmm31 := coeffn(mmm31,eps,0)$
nmm31 := coeffn(nmm31,eps,0)$
pmm31 := coeffn(pmm31,eps,0)$

lmm31_target := - m*(m-3)*(df(H3,r)-i*df(V3,r))$
mmm31_target := - (1/3)*m    *nmqm3*((1-sh)*(df(H3,r)-i*df(V3,r)) - (3*3-1)*(H3-i*V3)/r)$
nmm31_target := - (1/3)*(m-3)*nmq0 *((1-sh)*(df(H3,r)-i*df(V3,r)) - (3*3-1)*(H3-i*V3)/r)$
pmm31_target := - nmqm3*nmq0*(df(H3,r)-i*df(V3,r))$

write "Residuals of j=-3 coupling coefficients:";
dlmm31 := lmm31 - lmm31_target;
dmmm31 := mmm31 - mmm31_target;
dnmm31 := nmm31 - nmm31_target;
dpmm31 := pmm31 - pmm31_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=+4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmp41 := av(n*bmp41)$
mmp41 := av(cmp41 + lamp4*bmp41)$
nmp41 := av(dmp41 - lam00*bmp41)$
pmp41 := av(emp41 + lamp4*dmp41 - lam00*mmp41)/n$

lmp41 := coeffn(lmp41,eps,0)$
mmp41 := coeffn(mmp41,eps,0)$
nmp41 := coeffn(nmp41,eps,0)$
pmp41 := coeffn(pmp41,eps,0)$

lmp41_target := - m*(m+4)*(df(H4,r)+i*df(V4,r))$
mmp41_target :=   (1/4)*m    *nmqp4*((1-sh)*(df(H4,r)+i*df(V4,r)) - (4*4-1)*(H4+i*V4)/r)$
nmp41_target :=   (1/4)*(m+4)*nmq0 *((1-sh)*(df(H4,r)+i*df(V4,r)) - (4*4-1)*(H4+i*V4)/r)$
pmp41_target := - nmqp4*nmq0*(df(H4,r)+i*df(V4,r))$

write "Residuals of j=+4 coupling coefficients:";
dlmp41 := lmp41 - lmp41_target;
dmmp41 := mmp41 - mmp41_target;
dnmp41 := nmp41 - nmp41_target;
dpmp41 := pmp41 - pmp41_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate j=-4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmm41 := av(n*bmm41)$
mmm41 := av(cmm41 + lamm4*bmm41)$
nmm41 := av(dmm41 - lam00*bmm41)$
pmm41 := av(emm41 + lamm4*dmm41 - lam00*mmm41)/n$

lmm41 := coeffn(lmm41,eps,0)$
mmm41 := coeffn(mmm41,eps,0)$
nmm41 := coeffn(nmm41,eps,0)$
pmm41 := coeffn(pmm41,eps,0)$

lmm41_target := - m*(m-4)*(df(H4,r)-i*df(V4,r))$
mmm41_target := - (1/4)*m    *nmqm4*((1-sh)*(df(H4,r)-i*df(V4,r)) - (4*4-1)*(H4-i*V4)/r)$
nmm41_target := - (1/4)*(m-4)*nmq0 *((1-sh)*(df(H4,r)-i*df(V4,r)) - (4*4-1)*(H4-i*V4)/r)$
pmm41_target := - nmqm4*nmq0*(df(H4,r)-i*df(V4,r))$

write "Residuals of j=-4 coupling coefficients:";
dlmm41 := lmm41 - lmm41_target;
dmmm41 := mmm41 - mmm41_target;
dnmm41 := nmm41 - nmm41_target;
dpmm41 := pmm41 - pmm41_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate special j=0 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pm000 := av(sub(m=0,em0))/n$

pm000 := coeffn(pm000,eps,0)$
pm002 := coeffn(pm000,eps,2)$

res          := (r**2/q**2)*(2-sh)$
pm000_target := n**2*q**2 - q**2*df(res,r)/r - q**2*r*df(pp/r,r)$
pm002_target := 0$

write "Residuals of special j=0 coupling coefficients:";
dpm000 := pm000 - pm000_target;
dpm000 := pm002 - pm002_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate special case j=1 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmp10 := av(sub(m=+1,cmm11))$
mmm10 := av(sub(m=-1,cmp11))$
nm0p1 := av(sub(m=0, dmp11))$
nm0m1 := av(sub(m=0, dmm11))$
pm0p1 := av(sub(m=0, emp11) + sub(m=0,lamp1)*sub(m=0,dmp11))/n$
pm0m1 := av(sub(m=0, emm11) + sub(m=0,lamm1)*sub(m=0,dmm11))/n$

mmp10 := coeffn(mmp10,eps,0)$
mmm10 := coeffn(mmm10,eps,0)$
nm0p1 := coeffn(nm0p1,eps,0)$
nm0m1 := coeffn(nm0m1,eps,0)$
pm0p1 := coeffn(pm0p1,eps,0)$
pm0m1 := coeffn(pm0m1,eps,0)$

mmp10_target := sub(m=+1,mmm11_target) - (2-sh)*df(H1,r)$
mmm10_target := sub(m=-1,mmp11_target) + (2-sh)*df(H1,r)$
nm0p1_target := sub(m= 0,nmp11_target) + (2-sh)*df(H1,r)$
nm0m1_target := sub(m= 0,nmm11_target) - (2-sh)*df(H1,r)$
pm0p1_target := sub(m= 0,pmp11_target) - (2-sh)*(   n*q*pp*q**2 + (1-n*q)*(r + (1-sh)*df(H1,r)))$
pm0m1_target := sub(m= 0,pmm11_target) - (2-sh)*( - n*q*pp*q**2 + (1+n*q)*(r + (1-sh)*df(H1,r)))$

write "Residuals of special j=1 coupling coefficients:";
dmmp10 := mmp10 - mmp10_target;
dmmm10 := mmm10 - mmm10_target;
dnm0p1 := nm0p1 - nm0p1_target;
dnm0m1 := nm0m1 - nm0m1_target;
dpm0p1 := pm0p1 - pm0p1_target;
dpm0m1 := pm0m1 - pm0m1_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate special case j=2 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmp20 := av(sub(m=+2,cmm21))$
mmm20 := av(sub(m=-2,cmp21))$
nm0p2 := av(sub(m=0, dmp21))$
nm0m2 := av(sub(m=0, dmm21))$
pm0p2 := av(sub(m=0, emp21) + sub(m=0,lamp2)*sub(m=0,dmp21))/n$
pm0m2 := av(sub(m=0, emm21) + sub(m=0,lamm2)*sub(m=0,dmm21))/n$

mmp20 := coeffn(mmp20,eps,0)$
mmm20 := coeffn(mmm20,eps,0)$
nm0p2 := coeffn(nm0p2,eps,0)$
nm0m2 := coeffn(nm0m2,eps,0)$
pm0p2 := coeffn(pm0p2,eps,0)$
pm0m2 := coeffn(pm0m2,eps,0)$

mmp20_target := sub(m=+2,mmm21_target) - 2*(2-sh)*(df(H2,r)-i*df(V2,r))$
mmm20_target := sub(m=-2,mmp21_target) + 2*(2-sh)*(df(H2,r)+i*df(V2,r))$
nm0p2_target := sub(m= 0,nmp21_target) + 2*(2-sh)*(df(H2,r)+i*df(V2,r))$
nm0m2_target := sub(m= 0,nmm21_target) - 2*(2-sh)*(df(H2,r)-i*df(V2,r))$
pm0p2_target := sub(m= 0,pmp21_target) - (1/2)*(2-sh)*(2-n*q)*((1-sh)*(df(H2,r)+i*df(V2,r)) - (2*2-1)*(H2+i*V2)/r)$
pm0m2_target := sub(m= 0,pmm21_target) - (1/2)*(2-sh)*(2+n*q)*((1-sh)*(df(H2,r)-i*df(V2,r)) - (2*2-1)*(H2-i*V2)/r)$

write "Residuals of special j=2 coupling coefficients:";
dmmp20 := mmp20 - mmp20_target;
dmmm20 := mmm20 - mmm20_target;
dnm0p2 := nm0p2 - nm0p2_target;
dnm0m2 := nm0m2 - nm0m2_target;
dpm0p2 := pm0p2 - pm0p2_target;
dpm0m2 := pm0m2 - pm0m2_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate special case j=3 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmp30 := av(sub(m=+3,cmm31))$
mmm30 := av(sub(m=-3,cmp31))$
nm0p3 := av(sub(m=0, dmp31))$
nm0m3 := av(sub(m=0, dmm31))$
pm0p3 := av(sub(m=0, emp31) + sub(m=0,lamp3)*sub(m=0,dmp31))/n$
pm0m3 := av(sub(m=0, emm31) + sub(m=0,lamm3)*sub(m=0,dmm31))/n$

mmp30 := coeffn(mmp30,eps,0)$
mmm30 := coeffn(mmm30,eps,0)$
nm0p3 := coeffn(nm0p3,eps,0)$
nm0m3 := coeffn(nm0m3,eps,0)$
pm0p3 := coeffn(pm0p3,eps,0)$
pm0m3 := coeffn(pm0m3,eps,0)$

mmp30_target := sub(m=+3,mmm31_target) - 3*(2-sh)*(df(H3,r)-i*df(V3,r))$
mmm30_target := sub(m=-3,mmp31_target) + 3*(2-sh)*(df(H3,r)+i*df(V3,r))$
nm0p3_target := sub(m= 0,nmp31_target) + 3*(2-sh)*(df(H3,r)+i*df(V3,r))$
nm0m3_target := sub(m= 0,nmm31_target) - 3*(2-sh)*(df(H3,r)-i*df(V3,r))$
pm0p3_target := sub(m= 0,pmp31_target) - (1/3)*(2-sh)*(3-n*q)*((1-sh)*(df(H3,r)+i*df(V3,r)) - (3*3-1)*(H3+i*V3)/r)$
pm0m3_target := sub(m= 0,pmm31_target) - (1/3)*(2-sh)*(3+n*q)*((1-sh)*(df(H3,r)-i*df(V3,r)) - (3*3-1)*(H3-i*V3)/r)$

write "Residuals of special j=3 coupling coefficients:";
dmmp30 := mmp30 - mmp30_target;
dmmm30 := mmm30 - mmm30_target;
dnm0p3 := nm0p3 - nm0p3_target;
dnm0m3 := nm0m3 - nm0m3_target;
dpm0p3 := pm0p3 - pm0p3_target;
dpm0m3 := pm0m3 - pm0m3_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate special case j=4 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmp40 := av(sub(m=+4,cmm41))$
mmm40 := av(sub(m=-4,cmp41))$
nm0p4 := av(sub(m=0, dmp41))$
nm0m4 := av(sub(m=0, dmm41))$
pm0p4 := av(sub(m=0, emp41) + sub(m=0,lamp4)*sub(m=0,dmp41))/n$
pm0m4 := av(sub(m=0, emm41) + sub(m=0,lamm4)*sub(m=0,dmm41))/n$

mmp40 := coeffn(mmp40,eps,0)$
mmm40 := coeffn(mmm40,eps,0)$
nm0p4 := coeffn(nm0p4,eps,0)$
nm0m4 := coeffn(nm0m4,eps,0)$
pm0p4 := coeffn(pm0p4,eps,0)$
pm0m4 := coeffn(pm0m4,eps,0)$

mmp40_target := sub(m=+4,mmm41_target) - 4*(2-sh)*(df(H4,r)-i*df(V4,r))$
mmm40_target := sub(m=-4,mmp41_target) + 4*(2-sh)*(df(H4,r)+i*df(V4,r))$
nm0p4_target := sub(m= 0,nmp41_target) + 4*(2-sh)*(df(H4,r)+i*df(V4,r))$
nm0m4_target := sub(m= 0,nmm41_target) - 4*(2-sh)*(df(H4,r)-i*df(V4,r))$
pm0p4_target := sub(m= 0,pmp41_target) - (1/4)*(2-sh)*(4-n*q)*((1-sh)*(df(H4,r)+i*df(V4,r)) - (4*4-1)*(H4+i*V4)/r)$
pm0m4_target := sub(m= 0,pmm41_target) - (1/4)*(2-sh)*(4+n*q)*((1-sh)*(df(H4,r)-i*df(V4,r)) - (4*4-1)*(H4-i*V4)/r)$

write "Residuals of special j=4 coupling coefficients:";
dmmp40 := mmp40 - mmp40_target;
dmmm40 := mmm40 - mmm40_target;
dnm0p4 := nm0p4 - nm0p4_target;
dnm0m4 := nm0m4 - nm0m4_target;
dpm0p4 := pm0p4 - pm0p4_target;
dpm0m4 := pm0m4 - pm0m4_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate equilibrium quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g2p0 := coeffn(g2p,eps,0)$
g2p2 := coeffn(g2p,eps,2)$
g4p  := coeffn(g4p,eps,0)$

sp3 :=   df(H1,r)**2
       + df(H2,r)**2 + 2*(2*2-1)*df(H2,r)*H2/r - (2*2-1)*(H2/r)**2
       + df(H3,r)**2 + 2*(3*3-1)*df(H3,r)*H3/r - (3*3-1)*(H3/r)**2
       + df(H4,r)**2 + 2*(4*4-1)*df(H4,r)*H4/r - (4*4-1)*(H4/r)**2
       + df(V2,r)**2 + 2*(2*2-1)*df(V2,r)*V2/r - (2*2-1)*(V2/r)**2
       + df(V3,r)**2 + 2*(3*3-1)*df(V3,r)*V3/r - (3*3-1)*(V3/r)**2
       + df(V4,r)**2 + 2*(4*4-1)*df(V4,r)*V4/r - (4*4-1)*(V4/r)**2$

g2p_target  := - pp - r*(2-sh)/q**2$
g4p_target  := - (3*r*r/2 - 2*r*df(H1,r) + sp3)*r/q**2
               + (- g2 - 3*r**2/4 + r**2/q**2 + H1 + sp1)*(2-sh)*r/q**2
	       + pp * (g2 + r**2/2 + r**2/q**2 - 2*H1 - 3*r*df(H1,r))$
H1pp_target := - (3-2*sh)*df(H1,r)/r - 1 + 2*pp*q**2/r$
H2pp_target := - (3-2*sh)*df(H2,r)/r + (2*2-1)*H2/r**2$
H3pp_target := - (3-2*sh)*df(H3,r)/r + (3*3-1)*H3/r**2$
H4pp_target := - (3-2*sh)*df(H4,r)/r + (4*4-1)*H4/r**2$
V2pp_target := - (3-2*sh)*df(V2,r)/r + (2*2-1)*V2/r**2$
V3pp_target := - (3-2*sh)*df(V3,r)/r + (3*3-1)*V3/r**2$
V4pp_target := - (3-2*sh)*df(V4,r)/r + (4*4-1)*V4/r**2$

sp4 :=    df(H1,r)**2 
       +  df(H2,r)**2 + df(H3,r)**2 +  df(H4,r)**2
       +  df(V2,r)**2 + df(V3,r)**2 +  df(V4,r)**2$
sp5 :=    (1*1-1)*df(H1,r)*H1/r
       +  (2*2-1)*df(H2,r)*H2/r + (3*3-1)*df(H3,r)*H3/r +  (4*4-1)*df(H4,r)*H4/r
       +  (2*2-1)*df(V2,r)*V2/r + (3*3-1)*df(V3,r)*V3/r +  (4*4-1)*df(V4,r)*V4/r$
sp6 :=    (1*1-1)*H1**2/r**2
       +  (2*2-1)*H2**2/r**2 + (3*3-1)*H3**2/r**2 +  (4*4-1)*H4**2/r**2
       +  (2*2-1)*V2**2/r**2 + (3*3-1)*V3**2/r**2 +  (4*4-1)*V4**2/r**2$

kmp0_target := - (2 - sh)/m$
kmp1_target := 0$
kmp2_target := (r*pp - (3/2)*r**2 + 2*r*df(H1,r) - sp4 - 2*sp5 + sp6
   + (2-sh) * (- (3/4)*r**2 + r**2/q**2 + H1 + (1/2)*(3*sp4 - sp6)))/m
   + (n*r/m**2) * (-pp*q - r*(2-sh)*(m-n*q)/m/q)$

write "Equilibrium residuals:";
dg2p  := g2p  - g2p_target;
dg4p  := g4p  - g4p_target;
dH1pp := H1pp - H1pp_target;
dH2pp := H2pp - H2pp_target;
dV2pp := V2pp - V2pp_target;
dH3pp := H3pp - H3pp_target;
dV3pp := V3pp - V3pp_target;
dH4pp := H4pp - H4pp_target;
dV4pp := V4pp - V4pp_target;

dkmp0 := kmp0 - kmp0_target;
dkmp1 := kmp1 - kmp1_target;
dkmp2 := kmp2 - kmp2_target;

out T$

;bye;
