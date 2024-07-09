% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recreation of Jim Hastie and Richard Fitzpatrick's T7 calculation
% Extended for TJ calculation
% Extension for limied higher order calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out "Higher.out"$

% %%%%%%%%%%%%%
linelength 120$
% %%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express inverse of 3rd order polynomial ex (in var)
% as another 3rd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure inv(ex,var)$
   begin scalar c0,c1,c2,c3,i0,i1,i2,i3$
   c0 := coeffn(ex,var,0)$
   c1 := coeffn(ex,var,1)/c0$
   c2 := coeffn(ex,var,2)/c0$
   c3 := coeffn(ex,var,3)/c0$
   i0 := 1/c0$
   i1 := -c1/c0$
   i2 := (c1**2 - c2)/c0$
   i3 := (2*c1*c2 - c3 - c1*c1*c1)/c0
   return i0 + var*i1 + var**2*i2 + var**3*i3$
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
korder eps, i$
order  eps, i$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that all quantities expressed as polynomials in eps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
factor eps, cos, sin$
on  rat$
on  revpri$
off ratpri$

% %%%%%%%%%%%%%%%%%%%%%%%%
% Truncate at order eps**3
% %%%%%%%%%%%%%%%%%%%%%%%%
let eps**4 = 0$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fundamental equilibrium quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depend H1,r$
depend H2,r$
depend H3,r$
depend H4,r$
depend H5,r$
depend H6,r$
depend V2,r$
depend V3,r$
depend V4,r$
depend V5,r$
depend V6,r$
depend p2,r$
depend x,r,w$
depend z,r,w$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express x and z in terms of flux-surface label r, geometric
% angle w, and previously defined quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x := - r*cos(w) + eps*H1
       		+ eps*H2*cos(w) + eps*H3*cos(2*w) + eps*H4*cos(3*w) + eps*H5*cos(4*w) + eps*H6*cos(5*w)
		+ eps*V2*sin(w) + eps*V3*sin(2*w) + eps*V4*sin(3*w) + eps*V5*sin(4*w) + eps*V6*sin(5*w)
		+ eps**2*p2*cos(w) + eps**3*p3*cos(w)$
z :=   r*sin(w) + eps*H2*sin(w) + eps*H3*sin(2*w) + eps*H4*sin(3*w) + eps*H5*sin(4*w) + eps*H6*sin(5*w)
                - eps*V2*cos(w) - eps*V3*cos(2*w) - eps*V4*cos(3*w) - eps*V5*cos(4*w) - eps*V6*cos(5*w)
       		- eps**2*p2*sin(w) - eps**3*p3*sin(w)$

% %%%%%%%%%%%%%%%%%
% Evaluate Jacobian
% %%%%%%%%%%%%%%%%%
xr := df(x,r)$
xw := df(x,w)$
zr := df(z,r)$
zw := df(z,w)$
j  := xw*zr - xr*zw$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make transformation to straight field-line coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1  := 1 + eps*x$
ix1 := inv(x1,eps)$
k   := j*ix1$
ik  := int(k,w)$
ik0 := (sub(w=2*pi, ik) - sub(w=0,ik)) /(2*pi)$

write "Relabeling parameter:";
p2 :=
   r**3/8 - r*H1/2
   - (2-1)*H2**2/(2*r) - (3-1)*H3**2/(2*r) - (4-1)*H4**2/(2*r) - (5-1)*H5**2/(2*r) - (6-1)*H6**2/(2*r)
   - (2-1)*V2**2/(2*r) - (3-1)*V3**2/(2*r) - (4-1)*V4**2/(2*r) - (5-1)*V5**2/(2*r) - (6-1)*V6**2/(2*r)$
p3 :=
   - df(H2,r)*r*r*r/4 + H2*r*r/4
   -   (df(H2,r)*H3 + df(H3,r)*H2)*r/2
   - 2*(df(H3,r)*H4 + df(H4,r)*H3)*r/2
   - 3*(df(H4,r)*H5 + df(H5,r)*H4)*r/2
   - 4*(df(H5,r)*H6 + df(H6,r)*H5)*r/2
   -   H2*H3/2
   - 2*H3*H4/2
   - 3*H4*H5/2
   - 4*H5*H6/2
   -   (df(V2,r)*V3 + df(V3,r)*V2)*r/2
   - 2*(df(V3,r)*V4 + df(V4,r)*V3)*r/2
   - 3*(df(V4,r)*V5 + df(V5,r)*V4)*r/2
   - 4*(df(V5,r)*V6 + df(V6,r)*V5)*r/2
   -   V2*V3/2
   - 2*V3*V4/2
   - 3*V4*V5/2
   - 4*V5*V6/2$

res2 := coeffn(ik0,eps,2);
res3 := coeffn(ik0,eps,3);

write "Straight field-line flux-surface label:";
rst := sqrt(2) * sqrt(int(ik0,r));

write "Secular variation residual:";
theta := ik/r$
rest := (sub(w=2*pi,theta) - sub(w=0,theta))/(2*pi) - 1;

write "theta := tt0 + eps*tt1 + eps**2*tt2";
theta := ik/r$
tt0   := coeffn(theta,eps,0)$
tt1   := coeffn(theta,eps,1)$
tt2   := coeffn(theta,eps,2)$

tt2_c1 := coeffn(tt2,cos(w),  1)$
tt2_c2 := coeffn(tt2,cos(2*w),1)$
tt2_c3 := coeffn(tt2,cos(3*w),1)$
tt2_c4 := coeffn(tt2,cos(4*w),1)$
tt2_c5 := coeffn(tt2,cos(5*w),1)$
tt2_c6 := coeffn(tt2,cos(6*w),1)$
tt2_c7 := coeffn(tt2,cos(7*w),1)$

tt2_s1 := coeffn(tt2,sin(w),  1)$
tt2_s2 := coeffn(tt2,sin(2*w),1)$
tt2_s3 := coeffn(tt2,sin(3*w),1)$
tt2_s4 := coeffn(tt2,sin(4*w),1)$
tt2_s5 := coeffn(tt2,sin(5*w),1)$
tt2_s6 := coeffn(tt2,sin(6*w),1)$
tt2_s7 := coeffn(tt2,sin(7*w),1)$

tt0_target := w$
tt1_target := r*sin(w)
   - (df(H1,r) - (1-1)*H1/r)*sin(w)  /1
   - (df(H2,r) - (2-1)*H2/r)*sin(2*w)/2
   - (df(H3,r) - (3-1)*H3/r)*sin(3*w)/3
   - (df(H4,r) - (4-1)*H4/r)*sin(4*w)/4
   - (df(H5,r) - (5-1)*H5/r)*sin(5*w)/5
   - (df(H6,r) - (6-1)*H6/r)*sin(6*w)/6
   + (df(V2,r) - (2-1)*V2/r)*cos(2*w)/2
   + (df(V3,r) - (3-1)*V3/r)*cos(3*w)/3
   + (df(V4,r) - (4-1)*V4/r)*cos(4*w)/4
   + (df(V5,r) - (5-1)*V5/r)*cos(5*w)/5
   + (df(V6,r) - (6-1)*V6/r)*cos(6*w)/6$

tt2_c1_target :=
   + V2/2
   + df(V2,r)*r/2
   + df(H1,r)*V2/r$
tt2_c2_target :=
   + df(V3,r)*r/4
   + df(H1,r)*V3/r$
tt2_c3_target :=
   + (- (3-2)*V2 + df(V2,r)*r)/(2*3)
   + (- (3-2)*V4 + df(V4,r)*r)/(2*3)
   + df(H1,r)*V4/r$
tt2_c4_target :=
   + (- (4-2)*V3 + df(V3,r)*r)/(2*4)
   + (- (4-2)*V5 + df(V5,r)*r)/(2*4)
   + df(H1,r)*V5/r$
tt2_c5_target :=
   + (- (5-2)*V4 + df(V4,r)*r)/(2*5)
   + (- (5-2)*V6 + df(V6,r)*r)/(2*5)
   + df(H1,r)*V6/r$
tt2_c6_target :=
   + (- (6-2)*V5 + df(V5,r)*r)/(2*6)$
tt2_c7_target :=
   + (- (7-2)*V6 + df(V6,r)*r)/(2*7)$

tt2_s1_target :=
   - H2/2
   - df(H2,r)*r/2
   - df(H1,r)*H2/r
   - 2*(df(H2,r)*H3 + df(V2,r)*V3)/(1*r)
   - 3*(df(H3,r)*H4 + df(V3,r)*V4)/(1*r)
   - 4*(df(H4,r)*H5 + df(V4,r)*V5)/(1*r)
   - 5*(df(H5,r)*H6 + df(V5,r)*V6)/(1*r)
   - 1*(df(H3,r)*H2 + df(V3,r)*V2)/(1*r)
   - 2*(df(H4,r)*H3 + df(V4,r)*V3)/(1*r)
   - 3*(df(H5,r)*H4 + df(V5,r)*V4)/(1*r)
   - 4*(df(H6,r)*H5 + df(V6,r)*V5)/(1*r)$
tt2_s2_target :=
   + r*r/4
   - df(H1,r)*r/4
   - df(H3,r)*r/4
   - df(H1,r)*H3/r
   - 3*(df(H2,r)*H4 + df(V2,r)*V4)/(2*r)
   - 4*(df(H3,r)*H5 + df(V3,r)*V5)/(2*r)
   - 5*(df(H4,r)*H6 + df(V4,r)*V6)/(2*r)
   - 1*(df(H4,r)*H2 + df(V4,r)*V2)/(2*r)
   - 2*(df(H5,r)*H3 + df(V5,r)*V3)/(2*r)
   - 3*(df(H6,r)*H4 + df(V6,r)*V4)/(2*r)$
tt2_s3_target :=
   + ((3-2)*H2 - df(H2,r)*r)/(2*3)
   + ((3-2)*H4 - df(H4,r)*r)/(2*3)
   - df(H1,r)*H4/r
   - 4*(df(H2,r)*H5 + df(V2,r)*V5)/(3*r)
   - 5*(df(H3,r)*H6 + df(V3,r)*V6)/(3*r)
   - 1*(df(H5,r)*H2 + df(V5,r)*V2)/(3*r)
   - 2*(df(H6,r)*H3 + df(V6,r)*V3)/(3*r)$
tt2_s4_target :=
   + ((4-2)*H3 - df(H3,r)*r)/(2*4)
   + ((4-2)*H5 - df(H5,r)*r)/(2*4)
   - df(H1,r)*H5/r
   - 5*(df(H2,r)*H6 + df(V2,r)*V6)/(4*r)
   - 1*(df(H6,r)*H2 + df(V6,r)*V2)/(4*r)$
tt2_s5_target :=
   + ((5-2)*H4 - df(H4,r)*r)/(2*5)
   + ((5-2)*H6 - df(H6,r)*r)/(2*5)
   - df(H1,r)*H6/r$
tt2_s6_target :=
   + ((6-2)*H5 - df(H5,r)*r)/(2*6)$
tt2_s7_target :=
   + ((7-2)*H6 - df(H6,r)*r)/(2*7)$

tt0_residual := tt0 - tt0_target;
tt1_residual := tt1 - tt1_target;

tt2_c1_residual := tt2_c1 - tt2_c1_target;
tt2_c2_residual := tt2_c2 - tt2_c2_target;
tt2_c3_residual := tt2_c3 - tt2_c3_target;
tt2_c4_residual := tt2_c4 - tt2_c4_target;
tt2_c5_residual := tt2_c5 - tt2_c5_target;
tt2_c6_residual := tt2_c6 - tt2_c6_target;
tt2_c7_residual := tt2_c7 - tt2_c7_target;

tt2_s1_residual := tt2_s1 - tt2_s1_target;
tt2_s2_residual := tt2_s2 - tt2_s2_target;
tt2_s3_residual := tt2_s3 - tt2_s3_target;
tt2_s4_residual := tt2_s4 - tt2_s4_target;
tt2_s5_residual := tt2_s5 - tt2_s5_target;
tt2_s6_residual := tt2_s6 - tt2_s6_target;
tt2_s7_residual := tt2_s7 - tt2_s7_target;

out T$

;bye;
