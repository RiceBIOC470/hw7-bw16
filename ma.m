function [ timer,xsol ] = ma( a,xo )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rhs2=@(t,x) a*x.*(1-x);
interval=[0 10]; 
sol=ode23(rhs2,interval,xo); 

xsol=sol.x;
ysol=sol.y;
plot(xsol,ysol);

timen=.99<ysol;
timep=min(ysol(timen));
forx=min(ind2sub(size(ysol),find(timep==ysol)));
timer=xsol(forx);


end
