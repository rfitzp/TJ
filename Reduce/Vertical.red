% ############################################################################
% Script to calculate coupling coefficients for vertical stability calculation
% ############################################################################

out "Vertical.out"$

% %%%%%%%%%%%%%
linelength 120$
% %%%%%%%%%%%%%

depend H1,r$
depend S1,r$
depend S2,r$
depend g2,r$
depend p2p,r$
depend q,r$
depend s,r$

df(q,r)    := q * s /r$
df(g2,r)   := - p2p - r * (2 - s) /q/q$
g4p        := g2 * (p2p - r * (2  - s) /q/q) - (r/q) * S2 + p2p * (r*r/2 + r*r/q/q - 2*H1 - 3*r*df(H1,r))$
df(H1,r,2) := - (3 - 2*s) * df(H1,r) /r - 1 + 2*p2p*q*q/r$

alpp0 := p2p * q*q /r$
alpp2 := - 2 * alpp0 * g2$
alpg0 := q * df(g2,r) /r$
alpg2 := (q/r) * (- g2 * df(g2,r) + g4p)$
alpf0 := - s$
alpf2 := r * df(g2,r)$
J     := (2 - s) /q$

% %%%%%%%%%%%%%%
% Calculate dmm1
% %%%%%%%%%%%%%%
dmm1        := alpf0 * alpp0 + r * df(alpp0,r)$
dmm1_target := q*q * (df(r*p2p,r) - (2 - s) * p2p) /r$
dmm1_res    := dmm1 - dmm1_target;

% %%%%%%%%%%%%%%
% Calculate dmm0
% %%%%%%%%%%%%%%
dmm0        := - alpf0 * alpp0 - r * df(alpp0,r) - q * r * df(alpg0,r) + m*m$
dmm0_target := m*m + q * r * df(J,r)$
dmm0_res    := dmm0 - dmm0_target;

% %%%%%%%%%%%%%%
% Calculate dmm2
% %%%%%%%%%%%%%%
amm2 := - r*r/2 + r*df(H1,r) + 2*H1$
bmm2 := 7*r*r/4 - H1 - 3*r*df(H1,r) + S1$

dmm2 := - (alpf0*alpp0 + r*df(alpp0,r)) * amm2 - alpf0 * alpp2 - alpf2 * alpp0 - r * df(alpp2,r) - q*r * df(alpg2,r) - r*r *alpg0*alpg0 + m*m * bmm2$

dmm2_target := m*m * bmm2 - r*r*(2-s)*(2-s) /q/q + q*r*df(S2,r) - r*df(r*p2p,r) - 2*(1-s)*r*p2p
    + 2*r*p2p*q*q*(-2 + 3*p2p*q*q/r) + 2*df(H1,r)*q*q * (df(r*p2p,r) - 4*(1-s)*p2p)$

dmm2_res := dmm2 - dmm2_target;

out T$

;bye;
