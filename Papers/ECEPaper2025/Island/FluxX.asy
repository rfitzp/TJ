import graph;
import gsl;

size(750,500,IgnoreAspect);

real delta = 0.5;
real Omega;
real Y;

real rs (real z)
{
	return delta*cos(z*pi) /sqrt(8.);
}

real Xp (real z)
{
	real d2   = delta*delta;
	real cosx = cos (z*pi - d2*sin(z*pi));
	real Y2   = (Omega - cosx) /8.;

	if (Y2 > 0.)
   	   Y = sqrt(Y2);
	else
	   Y = 0.;

        return Y + delta*cos(z*pi) /sqrt(8.);
}

real Xm (real z)
{
	real d2   = delta*delta;
	real cosx = cos (z*pi - d2*sin(z*pi));
	real Y2   = (Omega - cosx)/8.;

	if (Y2 > 0.)
  	    Y = sqrt(Y2);
	else
	    Y = 0.;

	return - Y + delta*cos(z*pi) /sqrt(8.);
}

real X (real z)
{
	return Y + delta*cos(z*pi) /sqrt(8.);
}

pen s = black + solid + 1;

for (int i = 0; i < 20; ++i)
{
	if (i == 0)
	   Omega = -0.99;
	else
           Omega = 0.2*i - 1.;

        draw(graph(Xp,0.,4.,1000),s);
        draw(graph(Xm,0.,4.,1000),s);
}

s = solid + 2;
draw(graph(rs,0.,4.,1000),s);

s = solid + 2;
Omega = 1.;
draw(graph(Xp,0.,4.,1000),s);
draw(graph(Xm,0.,4.,1000),s);

s = dotted + 1;

for (int i = 0; i <= 20; ++i)
{
	Y = (-1. + 0.1*i);
	draw(graph(X,0.,4.,1000),s);
}

limits((0.,-1.), (4.,1.), Crop);     

for (int i = 0; i <= 40; ++i)
{
	real x = i*0.1*pi;
	real z = x;

	for (int j = 1; j < 40; ++j)
	{
		real n = j;
		z      = z + 2.*Jn(j,n*delta*delta)*sin(n*x)/n;
	}

	xequals (z/pi,s);
}

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$\zeta / \pi$",BottomTop,LeftTicks);
yaxis("$X$",LeftRight,RightTicks);
