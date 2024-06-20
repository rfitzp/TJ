% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recreation of Richard Fitzpatrick's small-r T7 calculation
% Extended for TJ calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out "Smallr.out"$

% %%%%%%%%%%%%%
linelength 120$
order eps$
korder r$
% %%%%%%%%%%%%%

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to return complex conjugate of expression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure cconj(ex)$
begin scalar c0,c1$
   c0 := coeffn(ex,i,0)$
   c1 := coeffn(ex,i,1)$
   return c0 - i*c1$
end$   

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expand functions of r about r=0
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depend q,r$
depend p,r$
depend H1,r$
depend H2,r$
depend V2,r$
depend H3,r$
depend V3,r$
depend H4,r$
depend V4,r$

q  := q0 + q2*r**2$
p  := p0 + p2*r**2$
pp := df(p,r)$
H1 := (1/8)*(4*p2*q0**2 - 1)*r**2$
H2 := H20*(1 + (q2/q0)*((2-1)/(2+1))*r**2)*r$
V2 := V20*(1 + (q2/q0)*((2-1)/(2+1))*r**2)*r$
H3 := H30*(1 + (q2/q0)*((3-1)/(3+1))*r**2)*r**2$
V3 := V30*(1 + (q2/q0)*((3-1)/(3+1))*r**2)*r**2$
H4 := H40*(1 + (q2/q0)*((4-1)/(4+1))*r**2)*r**3$
V4 := V40*(1 + (q2/q0)*((4-1)/(4+1))*r**2)*r**3$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate m,m coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmq0 := m - n*q$
iq   := inv(q,r)$
sh   := r*df(q,r)*iq$
r2q  := r**2*iq$
jcy  := df(r2q,r)/r$

lmm := m*m$
mmm := 0$
nmm := 0$
pmm := nmq0**2 + (nmq0/m)*(q*r)*df(jcy,r)$

lmm0 := coeffn(lmm,r,0)$
lmm2 := coeffn(lmm,r,2)$
mmm0 := coeffn(mmm,r,0)$
mmm2 := coeffn(mmm,r,2)$
nmm0 := coeffn(nmm,r,0)$
nmm2 := coeffn(nmm,r,2)$
pmm0 := coeffn(pmm,r,0)$
pmm2 := coeffn(pmm,r,2)$

on factor$
write "(m, m) coupling coefficients:";
lmm0 := lmm0;
lmm2 := lmm2;
mmm0 := mmm0;
mmm2 := mmm2;
nmm0 := nmm0;
nmm2 := nmm2;
pmm0 := pmm0;
pmm2 := pmm2;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,0 coupling coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iq2 := iq*iq$
res := r**2*iq2*(2-sh)$

l00 := 0;
m00 := 0;
n00 := 0;
p00 := n**2*q**2 - q**2*df(res,r)/r - q**2*r*df(pp/r,r)$

l000 := coeffn(l00,r,0)$
l002 := coeffn(l00,r,2)$
m000 := coeffn(m00,r,0)$
m002 := coeffn(m00,r,2)$
n000 := coeffn(n00,r,0)$
n002 := coeffn(n00,r,2)$
p000 := coeffn(p00,r,0)$
p002 := coeffn(p00,r,2)$

on factor$
write "(0, 0) coupling coefficients:";
l000 := l000;
l002 := l002;
m000 := m000;
m002 := m002;
n000 := n000;
n002 := n002;
p000 := p000;
p002 := p002;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate m,m+1 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmqp1 := m + 1 - n*q$

lmmp1 := - m*(m+1)*df(H1,r)$
mmmp1 := - m*nmq0*pp*q**2      + m*nmqp1   *(r + df(H1,r)*(1-sh))$
nmmp1 := - (m+1)*nmqp1*pp*q**2 + (m+1)*nmq0*(r + df(H1,r)*(1-sh))$
pmmp1 := - (1+sh)*pp*q**2      + nmq0*nmqp1*(r - df(H1,r))$

fac    := r$
lmmp10 := coeffn(lmmp1*fac,r,0)$
lmmp12 := coeffn(lmmp1*fac,r,2)$
mmmp10 := coeffn(mmmp1*fac,r,0)$
mmmp12 := coeffn(mmmp1*fac,r,2)$
nmmp10 := coeffn(nmmp1*fac,r,0)$
nmmp12 := coeffn(nmmp1*fac,r,2)$
pmmp10 := coeffn(pmmp1*fac,r,0)$
pmmp12 := coeffn(pmmp1*fac,r,2)$

on factor$
write "(m, m+1) coupling coefficients:";
lmmp10 := lmmp10;
lmmp12 := lmmp12;
mmnp10 := mmmp10;
mmmp12 := mmmp12;
nmmp10 := nmmp10;
nmmp12 := nmmp12;
pmmp10 := pmmp10;
pmmp12 := pmmp12;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,1 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l01 := sub(m=0,lmmp1)$
m01 := sub(m=0,mmmp1)$
n01 := sub(m=0,nmmp1) + (2-sh)*df(H1,r)$
p01 := sub(m=0,pmmp1) - (2-sh)*(n*q**3*pp + (1-n*q)*(r + (1-sh)*df(H1,r)))$

fac  := r$
l010 := coeffn(l01*fac,r,0)$
l012 := coeffn(l01*fac,r,2)$
m010 := coeffn(m01*fac,r,0)$
m012 := coeffn(m01*fac,r,2)$
n010 := coeffn(n01*fac,r,0)$
n012 := coeffn(n01*fac,r,2)$
p010 := coeffn(p01*fac,r,0)$
p012 := coeffn(p01*fac,r,2)$

on factor$
write "(0, 1) coupling coefficients:";
l010 := l010;
l012 := l012;
m010 := m010;
m012 := m012;
n010 := n010;
n012 := n012;
p010 := p010;
p012 := p012;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate -m,-m-1 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mx     := - m;
nmxq0  := mx - n*q$
nmxqm1 := mx - 1 - n*q$

lmxmm1 := - mx*(mx-1)*df(H1,r)$
mmxmm1 :=   mx*nmxq0*pp*q**2      - mx*nmxqm1   *(r + df(H1,r)*(1-sh))$
nmxmm1 :=   (mx-1)*nmxqm1*pp*q**2 - (mx-1)*nmxq0*(r + df(H1,r)*(1-sh))$
pmxmm1 := - (1+sh)*pp*q**2        + nmxq0*nmxqm1*(r - df(H1,r))$

fac     := r$
lmxmm10 := coeffn(lmxmm1*fac,r,0)$
lmxmm12 := coeffn(lmxmm1*fac,r,2)$
mmxmm10 := coeffn(mmxmm1*fac,r,0)$
mmxmm12 := coeffn(mmxmm1*fac,r,2)$
nmxmm10 := coeffn(nmxmm1*fac,r,0)$
nmxmm12 := coeffn(nmxmm1*fac,r,2)$
pmxmm10 := coeffn(pmxmm1*fac,r,0)$
pmxmm12 := coeffn(pmxmm1*fac,r,2)$

on factor$
write "(-m, -m-1) coupling coefficients:";
lmxmm10 := lmxmm10;
lmxmm12 := lmxmm12;
mmxnm10 := mmxmm10;
mmxmm12 := mmxmm12;
nmxmm10 := nmxmm10;
nmxmm12 := nmxmm12;
pmxmm10 := pmxmm10;
pmxmm12 := pmxmm12;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,-1 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l0m1 := sub(m=0,lmxmm1)$
m0m1 := sub(m=0,mmxmm1)$
n0m1 := sub(m=0,nmxmm1) - (2-sh)*df(H1,r)$
p0m1 := sub(m=0,pmxmm1) - (2-sh)*(- n*q**3*pp + (1+n*q)*(r + (1-sh)*df(H1,r)))$

fac   := r$
l0m10 := coeffn(l0m1*fac,r,0)$
l0m12 := coeffn(l0m1*fac,r,2)$
m0m10 := coeffn(m0m1*fac,r,0)$
m0m12 := coeffn(m0m1*fac,r,2)$
n0m10 := coeffn(n0m1*fac,r,0)$
n0m12 := coeffn(n0m1*fac,r,2)$
p0m10 := coeffn(p0m1*fac,r,0)$
p0m12 := coeffn(p0m1*fac,r,2)$

on factor$
write "(0, -1) coupling coefficients:";
l0m10 := l0m10;
l0m12 := l0m12;
m0m10 := m0m10;
m0m12 := m0m12;
n0m10 := n0m10;
n0m12 := n0m12;
p0m10 := p0m10;
p0m12 := p0m12;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate m,m+2 coupling coefficients normalized to 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmqp2 := m + 2 - n*q$

lmmp2 := - m*(m+2)*(df(H2,r)+i*df(V2,r))$
mmmp2 :=   (1/2)*m    *nmqp2*((1-sh)*(df(H2,r)+i*df(V2,r)) - (2*2-1)*(H2+i*V2)/r)$
nmmp2 :=   (1/2)*(m+2)*nmq0 *((1-sh)*(df(H2,r)+i*df(V2,r)) - (2*2-1)*(H2+i*V2)/r)$
pmmp2 := - nmqp2*nmq0*(df(H2,r)+i*df(V2,r))$

fac    := r**0$
lmmp20 := coeffn(lmmp2*fac,r,0)$
lmmp22 := coeffn(lmmp2*fac,r,2)$
mmmp20 := coeffn(mmmp2*fac,r,0)$
mmmp22 := coeffn(mmmp2*fac,r,2)$
nmmp20 := coeffn(nmmp2*fac,r,0)$
nmmp22 := coeffn(nmmp2*fac,r,2)$
pmmp20 := coeffn(pmmp2*fac,r,0)$
pmmp22 := coeffn(pmmp2*fac,r,2)$

on factor$
write "(m, m+2) coupling coefficients:";
lmmp20 := lmmp20;
lmmp22 := lmmp22;
mmnp20 := mmmp20;
mmmp22 := mmmp22;
nmmp20 := nmmp20;
nmmp22 := nmmp22;
pmmp20 := pmmp20;
pmmp22 := pmmp22;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,2 coupling coefficients normalized to 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l02 := sub(m=0,lmmp2)$
m02 := sub(m=0,mmmp2)$
n02 := sub(m=0,nmmp2) + 2*(2-sh)*(df(H2,r)+i*df(V2,r))$
p02 := sub(m=0,pmmp2) - (1/2)*(2-sh)*(2-n*q)*((1-sh)*(df(H2,r)+i*df(V2,r)) - (2*2-1)*(H2+i*V2)/r)$

fac  := r**0$
l020 := coeffn(l02*fac,r,0)$
l022 := coeffn(l02*fac,r,2)$
m020 := coeffn(m02*fac,r,0)$
m022 := coeffn(m02*fac,r,2)$
n020 := coeffn(n02*fac,r,0)$
n022 := coeffn(n02*fac,r,2)$
p020 := coeffn(p02*fac,r,0)$
p022 := coeffn(p02*fac,r,2)$

on factor$
write "(0, 2) coupling coefficients:";
l020 := l020;
l022 := l022;
m020 := m020;
m022 := m022;
n020 := n020;
n022 := n022;
p020 := p020;
p022 := p022;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate -m,-m-2 coupling coefficients normalized to 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmxqm2 := mx - 2 - n*q$

lmxmm2 := - mx*(mx-2)*(df(H2,r)-i*df(V2,r))$
mmxmm2 := - (1/2)*mx    *nmxqm2*((1-sh)*(df(H2,r)-i*df(V2,r)) - (2*2-1)*(H2-i*V2)/r)$
nmxmm2 := - (1/2)*(mx-2)*nmxq0 *((1-sh)*(df(H2,r)-i*df(V2,r)) - (2*2-1)*(H2-i*V2)/r)$
pmxmm2 := - nmxqm2*nmxq0*(df(H2,r)-i*df(V2,r))$

fac     := r**0$
lmxmm20 := coeffn(lmxmm2*fac,r,0)$
lmxmm22 := coeffn(lmxmm2*fac,r,2)$
mmxmm20 := coeffn(mmxmm2*fac,r,0)$
mmxmm22 := coeffn(mmxmm2*fac,r,2)$
nmxmm20 := coeffn(nmxmm2*fac,r,0)$
nmxmm22 := coeffn(nmxmm2*fac,r,2)$
pmxmm20 := coeffn(pmxmm2*fac,r,0)$
pmxmm22 := coeffn(pmxmm2*fac,r,2)$

on factor$
write "(-m, -m-2) coupling coefficients:";
lmxmm20 := lmxmm20;
lmxmm22 := lmxmm22;
mmxnm20 := mmxmm20;
mmxmm22 := mmxmm22;
nmxmm20 := nmxmm20;
nmxmm22 := nmxmm22;
pmxmm20 := pmxmm20;
pmxmm22 := pmxmm22;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,-2 coupling coefficients normalized to 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l0m2 := sub(m=0,lmxmm2)$
m0m2 := sub(m=0,mmxmm2)$
n0m2 := sub(m=0,nmxmm2) - 2*(2-sh)*(df(H2,r)-i*df(V2,r))$
p0m2 := sub(m=0,pmxmm2) - (1/2)*(2-sh)*(2+n*q)*((1-sh)*(df(H2,r)-i*df(V2,r)) - (2*2-1)*(H2-i*V2)/r)$

fac   := r**0$
l0m20 := coeffn(l0m2*fac,r,0)$
l0m22 := coeffn(l0m2*fac,r,2)$
m0m20 := coeffn(m0m2*fac,r,0)$
m0m22 := coeffn(m0m2*fac,r,2)$
n0m20 := coeffn(n0m2*fac,r,0)$
n0m22 := coeffn(n0m2*fac,r,2)$
p0m20 := coeffn(p0m2*fac,r,0)$
p0m22 := coeffn(p0m2*fac,r,2)$

on factor$
write "(0, -2) coupling coefficients:";
l0m20 := l0m20;
l0m22 := l0m22;
m0m20 := m0m20;
m0m22 := m0m22;
n0m20 := n0m20;
n0m22 := n0m22;
p0m20 := p0m20;
p0m22 := p0m22;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate m,m+3 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmqp3 := m + 3 - n*q$

lmmp3 := - m*(m+3)*(df(H3,r)+i*df(V3,r))$
mmmp3 :=   (1/3)*m    *nmqp3*((1-sh)*(df(H3,r)+i*df(V3,r)) - (3*3-1)*(H3+i*V3)/r)$
nmmp3 :=   (1/3)*(m+3)*nmq0 *((1-sh)*(df(H3,r)+i*df(V3,r)) - (3*3-1)*(H3+i*V3)/r)$
pmmp3 := - nmqp3*nmq0*(df(H3,r)+i*df(V3,r))$

fac    := r$
lmmp30 := coeffn(lmmp3*fac,r,0)$
lmmp32 := coeffn(lmmp3*fac,r,2)$
mmmp30 := coeffn(mmmp3*fac,r,0)$
mmmp32 := coeffn(mmmp3*fac,r,2)$
nmmp30 := coeffn(nmmp3*fac,r,0)$
nmmp32 := coeffn(nmmp3*fac,r,2)$
pmmp30 := coeffn(pmmp3*fac,r,0)$
pmmp32 := coeffn(pmmp3*fac,r,2)$

on factor$
write "(m, m+3) coupling coefficients:";
lmmp30 := lmmp30;
lmmp32 := lmmp32;
mmnp30 := mmmp30;
mmmp32 := mmmp32;
nmmp30 := nmmp30;
nmmp32 := nmmp32;
pmmp30 := pmmp30;
pmmp32 := pmmp32;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,3 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l03 := sub(m=0,lmmp3)$
m03 := sub(m=0,mmmp3)$
n03 := sub(m=0,nmmp3) + 3*(2-sh)*(df(H3,r)+i*df(V3,r))$
p03 := sub(m=0,pmmp3) - (1/3)*(2-sh)*(3-n*q)*((1-sh)*(df(H3,r)+i*df(V3,r)) - (3*3-1)*(H3+i*V3)/r)$

fac  := r$
l030 := coeffn(l03*fac,r,0)$
l032 := coeffn(l03*fac,r,2)$
m030 := coeffn(m03*fac,r,0)$
m032 := coeffn(m03*fac,r,2)$
n030 := coeffn(n03*fac,r,0)$
n032 := coeffn(n03*fac,r,2)$
p030 := coeffn(p03*fac,r,0)$
p032 := coeffn(p03*fac,r,2)$

on factor$
write "(0, 3) coupling coefficients:";
l030 := l030;
l032 := l032;
m030 := m030;
m032 := m032;
n030 := n030;
n032 := n032;
p030 := p030;
p032 := p032;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate -m,-m-3 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmxqm3 := mx - 3 - n*q$

lmxmm3 := - mx*(mx-3)*(df(H3,r)-i*df(V3,r))$
mmxmm3 := - (1/3)*mx    *nmxqm3*((1-sh)*(df(H3,r)-i*df(V3,r)) - (3*3-1)*(H3-i*V3)/r)$
nmxmm3 := - (1/3)*(mx-3)*nmxq0 *((1-sh)*(df(H3,r)-i*df(V3,r)) - (3*3-1)*(H3-i*V3)/r)$
pmxmm3 := - nmxqm3*nmxq0*(df(H3,r)-i*df(V3,r))$

fac     := r$
lmxmm30 := coeffn(lmxmm3*fac,r,0)$
lmxmm32 := coeffn(lmxmm3*fac,r,2)$
mmxmm30 := coeffn(mmxmm3*fac,r,0)$
mmxmm32 := coeffn(mmxmm3*fac,r,2)$
nmxmm30 := coeffn(nmxmm3*fac,r,0)$
nmxmm32 := coeffn(nmxmm3*fac,r,2)$
pmxmm30 := coeffn(pmxmm3*fac,r,0)$
pmxmm32 := coeffn(pmxmm3*fac,r,2)$

on factor$
write "(-m, -m-3) coupling coefficients:";
lmxmm30 := lmxmm30;
lmxmm32 := lmxmm32;
mmxnm30 := mmxmm30;
mmxmm32 := mmxmm32;
nmxmm30 := nmxmm30;
nmxmm32 := nmxmm32;
pmxmm30 := pmxmm30;
pmxmm32 := pmxmm32;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate 0,-3 coupling coefficients normalized to r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l0m3 := sub(m=0,lmxmm3)$
m0m3 := sub(m=0,mmxmm3)$
n0m3 := sub(m=0,nmxmm3) - 3*(2-sh)*(df(H3,r)-i*df(V3,r))$
p0m3 := sub(m=0,pmxmm3) - (1/3)*(2-sh)*(3+n*q)*((1-sh)*(df(H3,r)-i*df(V3,r)) - (3*3-1)*(H3-i*V3)/r)$

fac   := r$
l0m30 := coeffn(l0m3*fac,r,0)$
l0m32 := coeffn(l0m3*fac,r,2)$
m0m30 := coeffn(m0m3*fac,r,0)$
m0m32 := coeffn(m0m3*fac,r,2)$
n0m30 := coeffn(n0m3*fac,r,0)$
n0m32 := coeffn(n0m3*fac,r,2)$
p0m30 := coeffn(p0m3*fac,r,0)$
p0m32 := coeffn(p0m3*fac,r,2)$

on factor$
write "(0, -3) coupling coefficients:";
l0m30 := l0m30;
l0m32 := l0m32;
m0m30 := m0m30;
m0m32 := m0m32;
n0m30 := n0m30;
n0m32 := n0m32;
p0m30 := p0m30;
p0m32 := p0m32;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve mth set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%
psim  := (1   + am2*r**2)*r**m$
zm    := (bm0 + bm2*r**2)*r**m$

fac  := r**(-m)$
eq1m := fac*(r*df(psim,r) - lmm*inv(nmq0,r)*zm)$
eq2m := fac*(r*df(zm,r)   - (pmm*psim - n*q*sh*zm)*inv(nmq0,r))$

% Lowest order solution
eq1m0 := coeffn(eq1m,r,0)$
eq2m0 := coeffn(eq2m,r,0)$

c0 := coeffn(eq1m0,bm0,0)$
c1 := coeffn(eq1m0,bm0,1)$

on factor$
bm0 := -c0/c1$
write "mth solution:";
am0 := 1;
bm0 := bm0;
eq1m0;
eq2m0;
off factor$

% Higher order solution
eq1m2 := coeffn(eq1m,r,2)$
eq2m2 := coeffn(eq2m,r,2)$

a11 := coeffn(eq1m2,am2,1)$
a12 := coeffn(eq1m2,bm2,1)$
a21 := coeffn(eq2m2,am2,1)$
a22 := coeffn(eq2m2,bm2,1)$

b1 := a11*am2 + a12*bm2 - eq1m2$
b2 := a21*am2 + a22*bm2 - eq2m2$

deta := a11*a22-a12*a21$

on factor$
am2 := (a22*b1 - a12*b2)/deta;
bm2 := (a22*b2 - a21*b1)/deta;
eq1m2;
eq2m2;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (m+1)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psim1 := am1*r**(m+1)$
zm1   := bm1*r**(m+1)$

lm1m1 := sub(m=kk+1,lmm)$
lm1m1 := sub(kk=m,lm1m1)$
pm1m1 := sub(m=kk+1,pmm)$
pm1m1 := sub(kk=m,pm1m1)$

fac   := r**(-(m+1))$
eq1m1 := fac*(r*df(psim1,r) - lm1m1*inv(nmqp1,r)*zm1
   - eps * (cconj(lmmp1)*zm - cconj(nmmp1)*psim)*inv(nmq0,r))$
eq2m1:= fac*(r*df(zm1,r) - (pm1m1*psim1 - n*q*sh*zm1)*inv(nmqp1,r)
   - eps * (- cconj(mmmp1)*zm + cconj(pmmp1)*psim)*inv(nmq0,r))$

eq1m10 := coeffn(eq1m1,r,0)$
eq2m10 := coeffn(eq2m1,r,0)$

a11 := coeffn(eq1m10,am1,1)$
a12 := coeffn(eq1m10,bm1,1)$
a21 := coeffn(eq2m10,am1,1)$
a22 := coeffn(eq2m10,bm1,1)$

b1 := a11*am1 + a12*bm1 - eq1m10$
b2 := a21*am1 + a22*bm1 - eq2m10$

deta := a11*a22-a12*a21$

on factor$
write "(m+1)th solution:";
am1 := b1/a11;
bm1 := 0;
eq1m10;
eq2m10;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (m+2)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psim2 := am22*r**(m+2)$
zm2   := bm22*r**(m+2)$

lm2m1 := sub(m=kk+2,lmm)$
lm2m1 := sub(kk=m,lm2m1)$
pm2m1 := sub(m=kk+2,pmm)$
pm2m1 := sub(kk=m,pm2m1)$

fac   := r**(-(m+2))$
eq1m2 := fac*(r*df(psim2,r) - lm2m1*inv(nmqp2,r)*zm2
   - eps * (cconj(lmmp2)*zm - cconj(nmmp2)*psim)*inv(nmq0,r))$
eq2m2:= fac*(r*df(zm2,r) - (pm2m1*psim2 - n*q*sh*zm2)*inv(nmqp2,r)
   - eps * (- cconj(mmmp2)*zm + cconj(pmmp2)*psim)*inv(nmq0,r))$

eq1m20 := coeffn(eq1m2,r,0)$
eq2m20 := coeffn(eq2m2,r,0)$

a11 := coeffn(eq1m20,am22,1)$
a12 := coeffn(eq1m20,bm22,1)$
a21 := coeffn(eq2m20,am22,1)$
a22 := coeffn(eq2m20,bm22,1)$

b1 := a11*am22 + a12*bm22 - eq1m20$
b2 := a21*am22 + a22*bm22 - eq2m20$

on factor$
write "(m+2)th solution:";
am22 := b1/a11;
bm22 := 0;
eq1m20;
eq2m20;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (m+3)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psim3 := am33*r**(m+3)$
zm3   := bm33*r**(m+3)$

lm3m1 := sub(m=kk+3,lmm)$
lm3m1 := sub(kk=m,lm3m1)$
pm3m1 := sub(m=kk+3,pmm)$
pm3m1 := sub(kk=m,pm3m1)$

fac   := r**(-(m+3))$
eq1m3 := fac*(r*df(psim3,r) - lm3m1*inv(nmqp3,r)*zm3
   - eps * (cconj(lmmp3)*zm - cconj(nmmp3)*psim)*inv(nmq0,r))$
eq2m3:= fac*(r*df(zm3,r) - (pm3m1*psim3 - n*q*sh*zm3)*inv(nmqp3,r)
   - eps * (- cconj(mmmp3)*zm + cconj(pmmp3)*psim)*inv(nmq0,r))$

eq1m30 := coeffn(eq1m3,r,0)$
eq2m30 := coeffn(eq2m3,r,0)$

a11 := coeffn(eq1m30,am33,1)$
a12 := coeffn(eq1m30,bm33,1)$
a21 := coeffn(eq2m30,am33,1)$
a22 := coeffn(eq2m30,bm33,1)$

b1 := a11*am33 + a12*bm33 - eq1m30$
b2 := a21*am33 + a22*bm33 - eq2m30$

on factor$
write "(m+3)th solution:";
am33 := b1/a11;
bm33 := 0;
eq1m30;
eq2m30;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve 0th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%
psi0 := a00 + a02*r**2$
z0   := b00 + b02*r**2$

znq := -n*q$
fac := r**0$
eq10 := fac*(r*df(psi0,r) - l00*inv(znq,r)*z0)$
eq20 := fac*(r*df(z0,r) - (p00*psi0 - n*q*sh*z0)*inv(znq,r))$

% Zeroth order solution
eq100 := coeffn(eq10,r,0)$
eq200 := coeffn(eq20,r,0)$

c0 := coeffn(eq200,a00,0)$
c1 := coeffn(eq200,a00,1)$

on factor$
write "0th solution:";
a00 := -c0/c1;
b00 := 1;
eq100;
eq200;
off factor$

% Higher order solution
eq102 := coeffn(eq10,r,2)$
eq202 := coeffn(eq20,r,2)$

c0  := coeffn(eq102,a02,0)$
c1  := coeffn(eq102,a02,1)$
on factor$
a02 := -c0/c1;
off factor$
c0  := coeffn(eq202,b02,0)$
c1  := coeffn(eq202,b02,1)$
on factor$
b02 := -c0/c1;
eq102;
eq202;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (0+1)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi01 := a01*r**1$
z01   := b01*r**1$

l11   := sub(m=1,lmm)$
p11   := sub(m=1,pmm)$
onq   := sub(m=1,nmq0)$

fac   := r**(-1)$
eq101 := fac*(r*df(psi101,r) - l11*inv(onq,r)*z01
   - eps * (cconj(l01)*z0 - cconj(n01)*psi0)*inv(znq,r))$
eq201:= fac*(r*df(z01,r) - (p11*psi01 - n*q*sh*z01)*inv(onq,r)
   - eps * (- cconj(m01)*z0 + cconj(p01)*psi0)*inv(znq,r))$

eq1010 := coeffn(eq101,r,0)$
eq2010 := coeffn(eq201,r,0)$

a11 := coeffn(eq1010,a01,1)$
a21 := coeffn(eq1010,b01,1)$
a12 := coeffn(eq2010,a01,1)$
a22 := coeffn(eq2010,b01,1)$

b1 := a11*a01 + a12*b01 - eq1010$
b2 := a21*a01 + a22*b01 - eq2010$

on factor$
write "(0+1)th solution:";
a01 := 0;
b01 := 0;
eq1010;
eq2010;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (0+2)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi02 := a002*r**1$
z02   := b002*r**1$

l22   := sub(m=2,lmm)$
p22   := sub(m=2,pmm)$
onq   := sub(m=2,nmq0)$

fac   := r**(-2)$
eq101 := fac*(r*df(psi102,r) - l22*inv(onq,r)*z01
   - eps * (cconj(l02)*z0 - cconj(n02)*psi0)*inv(znq,r))$
eq201:= fac*(r*df(z02,r) - (p22*psi02 - n*q*sh*z02)*inv(onq,r)
   - eps * (- cconj(m02)*z0 + cconj(p02)*psi0)*inv(znq,r))$

eq1020 := coeffn(eq102,r,0)$
eq2020 := coeffn(eq202,r,0)$

on factor$
write "(0+2)th solution:";
a002 := 0;
b002:= 0;
eq1020;
eq2020;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve -mth set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
psimm  := (1    + amm2*r**2)*r**m$
zmm    := (bmm0 + bmm2*r**2)*r**m$

lmmmm := sub(m=kk,lmm)$
lmmmm := sub(kk=-m,lmmmm)$
pmmmm := sub(m=kk,pmm)$
pmmmm := sub(kk=-m,pmmmm)$

fac   := r**(-m)$
eq1mm := fac*(r*df(psimm,r) - lmmmm*inv(nmxq0,r)*zmm)$
eq2mm := fac*(r*df(zmm,r)   - (pmmmm*psimm - n*q*sh*zmm)*inv(nmxq0,r))$

% Lowest order solution
eq1mm0 := coeffn(eq1mm,r,0)$
eq2mm0 := coeffn(eq2mm,r,0)$

c0 := coeffn(eq1mm0,bmm0,0)$
c1 := coeffn(eq1mm0,bmm0,1)$

on factor$
bmm0 := -c0/c1$
write "-mth solution:";
amm0 := 1;
bmm0 := bmm0;
eq1mm0;
eq2mm0;
off factor$

% Higher order solution
eq1mm2 := coeffn(eq1mm,r,2)$
eq2mm2 := coeffn(eq2mm,r,2)$

a11 := coeffn(eq1mm2,amm2,1)$
a12 := coeffn(eq1mm2,bmm2,1)$
a21 := coeffn(eq2mm2,amm2,1)$
a22 := coeffn(eq2mm2,bmm2,1)$

b1 := a11*amm2 + a12*bmm2 - eq1mm2$
b2 := a21*amm2 + a22*bmm2 - eq2mm2$

deta := a11*a22-a12*a21$

on factor$
amm2 := (a22*b1 - a12*b2)/deta;
bmm2 := (a22*b2 - a21*b1)/deta;
eq1mm2;
eq2mm2;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (-m-1)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psimm11 := amm11*r**(m+1)$
zmm11   := bmm11*r**(m+1)$

lmm11mm11 := sub(m=kk-1,lmm)$
lmm11mm11 := sub(kk=-m,lmm11mm11)$
pmm11mm11 := sub(m=kk-1,pmm)$
pmm11mm11 := sub(kk=-m,pmm11mm11)$

fac     := r**(-(m+1))$
eq1mm11 := fac*(r*df(psimm11,r) - lmm11mm11*inv(nmxqm1,r)*zmm11
   - eps * (cconj(lmxmm1)*zmm - cconj(nmxmm1)*psimm)*inv(nmxq0,r))$
eq2mm11 := fac*(r*df(zmm11,r) - (pmm11mm11*psimm11 - n*q*sh*zmm11)*inv(nmxqm1,r)
   - eps * (- cconj(mmxmm1)*zmm + cconj(pmxmm1)*psimm)*inv(nmxq0,r))$

eq1mm110 := coeffn(eq1mm11,r,0)$
eq2mm110 := coeffn(eq2mm11,r,0)$

a11 := coeffn(eq1mm110,amm11,1)$
a12 := coeffn(eq1mm110,bmm11,1)$
a21 := coeffn(eq2mm110,amm11,1)$
a22 := coeffn(eq2mm110,bmm11,1)$

b1 := a11*amm11 + a12*bmm11 - eq1mm110$
b2 := a21*amm11 + a22*bmm11 - eq2mm110$

deta := a11*a22-a12*a21$

on factor$
	 write "(-m-1)th solution:";
amm11 := b1/a11;
bmm11 := 0;
eq1mm110;
eq2mm110;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (-m-2)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psimm21 := amm21*r**(m+2)$
zmm21   := bmm21*r**(m+2)$

lmm11mm21 := sub(m=kk-2,lmm)$
lmm11mm21 := sub(kk=-m,lmm11mm21)$
pmm11mm21 := sub(m=kk-2,pmm)$
pmm11mm21 := sub(kk=-m,pmm11mm21)$

fac     := r**(-(m+2))$
eq1mm21 := fac*(r*df(psimm21,r) - lmm11mm21*inv(nmxqm2,r)*zmm21
   - eps * (cconj(lmxmm2)*zmm - cconj(nmxmm2)*psimm)*inv(nmxq0,r))$
eq2mm21 := fac*(r*df(zmm21,r) - (pmm11mm21*psimm21 - n*q*sh*zmm21)*inv(nmxqm2,r)
   - eps * (- cconj(mmxmm2)*zmm + cconj(pmxmm2)*psimm)*inv(nmxq0,r))$

eq1mm120 := coeffn(eq1mm21,r,0)$
eq2mm120 := coeffn(eq2mm21,r,0)$

a11 := coeffn(eq1mm120,amm21,1)$
a12 := coeffn(eq1mm120,bmm21,1)$
a21 := coeffn(eq2mm120,amm21,1)$
a22 := coeffn(eq2mm120,bmm21,1)$

b1 := a11*amm21 + a12*bmm21 - eq1mm120$
b2 := a21*amm21 + a22*bmm21 - eq2mm120$

deta := a11*a22-a12*a21$

on factor$
	 write "(-m-2)th solution:";
amm21 := b1/a11;
bmm21 := 0;
eq1mm120;
eq2mm120;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (-m-3)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psimm31 := amm31*r**(m+3)$
zmm31   := bmm31*r**(m+3)$

lmm11mm31 := sub(m=kk-3,lmm)$
lmm11mm31 := sub(kk=-m,lmm11mm31)$
pmm11mm31 := sub(m=kk-3,pmm)$
pmm11mm31 := sub(kk=-m,pmm11mm31)$

fac     := r**(-(m+3))$
eq1mm31 := fac*(r*df(psimm31,r) - lmm11mm31*inv(nmxqm3,r)*zmm31
   - eps * (cconj(lmxmm3)*zmm - cconj(nmxmm3)*psimm)*inv(nmxq0,r))$
eq2mm31 := fac*(r*df(zmm31,r) - (pmm11mm31*psimm31 - n*q*sh*zmm31)*inv(nmxqm3,r)
   - eps * (- cconj(mmxmm3)*zmm + cconj(pmxmm3)*psimm)*inv(nmxq0,r))$

eq1mm130 := coeffn(eq1mm31,r,0)$
eq2mm130 := coeffn(eq2mm31,r,0)$

a11 := coeffn(eq1mm130,amm31,1)$
a12 := coeffn(eq1mm130,bmm31,1)$
a21 := coeffn(eq2mm130,amm31,1)$
a22 := coeffn(eq2mm130,bmm31,1)$

b1 := a11*amm31 + a12*bmm31 - eq1mm130$
b2 := a21*amm31 + a22*bmm31 - eq2mm130$

deta := a11*a22-a12*a21$

on factor$
write "(-m-3)th solution:";
amm31 := b1/a11;
bmm31 := 0;
eq1mm130;
eq2mm120;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (0-1)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi011 := a011*r**1$
z011   := b011*r**1$

l1111  := sub(m=-1,lmm)$
p1111  := sub(m=-1,pmm)$
oonq   := sub(m=-1,nmq0)$

fac    := r**(-1)$
eq1011 := fac*(r*df(psi011,r) - l1111*inv(oonq,r)*z011
   - eps * (cconj(l0m1)*z0 - cconj(n0m1)*psi0)*inv(znq,r))$
eq2011:= fac*(r*df(z01,r) - (p11*psi01 - n*q*sh*z01)*inv(onq,r)
   - eps * (- cconj(m0m1)*z0 + cconj(p0m1)*psi0)*inv(znq,r))$

eq10110 := coeffn(eq1011,r,0)$
eq20110 := coeffn(eq2011,r,0)$

a11 := coeffn(eq10110,a011,1)$
a21 := coeffn(eq10110,b011,1)$
a12 := coeffn(eq20110,a011,1)$
a22 := coeffn(eq20110,b011,1)$

b1 := a11*a011 + a12*b011 - eq10110$
b2 := a21*a011 + a22*b011 - eq20110$

on factor$
write "(0-1)th solution:";
a011 := 0;
b011 := 0;
eq10110;
eq20110;
off factor$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve (0-2)th set of equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi022 := a022*r**2$
z022   := b022*r**2$

l2222  := sub(m=-2,lmm)$
p2222  := sub(m=-2,pmm)$
oonq   := sub(m=-2,nmq0)$

fac    := r**(-2)$
eq1022 := fac*(r*df(psi022,r) - l2222*inv(oonq,r)*z022
   - eps * (cconj(l0m2)*z0 - cconj(n0m2)*psi0)*inv(znq,r))$
eq2022:= fac*(r*df(z02,r) - (p22*psi02 - n*q*sh*z02)*inv(onq,r)
   - eps * (- cconj(m0m2)*z0 + cconj(p0m2)*psi0)*inv(znq,r))$

eq10220 := coeffn(eq1022,r,0)$
eq20220 := coeffn(eq2022,r,0)$

a11 := coeffn(eq10220,a022,1)$
a21 := coeffn(eq10220,b022,1)$
a12 := coeffn(eq20220,a022,1)$
a22 := coeffn(eq20220,b022,1)$

b1 := a11*a022 + a12*b022 - eq10220$
b2 := a21*a022 + a22*b022 - eq20220$

on factor$
write "(0-2)th solution:";
a022 := 0;
b022 := 0;
eq10220;
eq20220;
off factor$

out T$
;bye;
