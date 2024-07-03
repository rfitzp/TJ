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
		+ eps**2*p2*cos(w) + eps**3*p3*cos(w)$
z :=   r*sin(w) + eps*H2*sin(w) + eps*H3*sin(2*w) + eps*H4*sin(3*w)
                - eps*V2*cos(w) - eps*V3*cos(2*w) - eps*V4*cos(3*w)
       		- eps**2*p2*sin(w) - eps**3*p3*sin(w)$

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
p3   :=
   - df(H2,r)*r*r*r/4 + H2*r*r/4
   - df(H2,r)*H3*r/2 - df(H3,r)*H2*r/2 - df(H3,r)*H4*r - df(H4,r)*H3*r
   - H2*H3/2         - H3*H4
   - df(V2,r)*V3*r/2 - df(V3,r)*V2*r/2 - df(V3,r)*V4*r - df(V4,r)*V3*r
   - V2*V3/2         - V3*V4$
res2 := coeffn(ik0,eps,2);
res3 := coeffn(ik0,eps,3);

write "Straight field-line flux-surface label:";
rst := sqrt(2) * sqrt(int(ik0,r));

write "theta := tt0 + eps*tt1 + eps**2*tt2";
theta := ik/r$
tt0   := coeffn(theta,eps,0)$
tt1   := coeffn(theta,eps,1)$
tt2   := coeffn(theta,eps,2)$

tt0_target := w$
tt1_target := r * sin(w)
   - (df(H1, r) - (1-1)*H1/r) * sin(1*w)/1
   - (df(H2, r) - (2-1)*H2/r) * sin(2*w)/2
   - (df(H3, r) - (3-1)*H3/r) * sin(3*w)/3
   - (df(H4, r) - (4-1)*H4/r) * sin(4*w)/4
   + (df(V2, r) - (2-1)*V2/r) * cos(2*w)/2
   + (df(V3, r) - (3-1)*V3/r) * cos(3*w)/3
   + (df(V4, r) - (4-1)*V4/r) * cos(4*w)/4$
tt2_target := 0$

tt0_residual := tt0 - tt0_target;
tt1_residual := tt1 - tt1_target;
tt2_residual := tt2 - tt2_target;

write "Secular variation residual:";
rest := (sub(w=2*pi,theta) - sub(w=0,theta))/(2*pi) - 1;

out T$

;bye;
