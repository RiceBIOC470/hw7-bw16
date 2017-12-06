%HW7
%GB comments
1a 90 mentioning that the stable point at 1 is the carrying capacity of the population. Saying it stops once it reaches the limit is not accurate. It slows down and/or will be corrected (if values exceed 1) as the population converges towards 1. 
1b. 40 There is no explanation in how A could be impacting the system. It does NOT affect the stable points but it does affect how quickly the system is corrected towards the stable points. 
1c.  100
1d. 40 Graph is not correct and there is no explanation of your results 
2a 0 This is wrong. This is correct [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
2b. 30 so many issues with your line of code. you don’t define A and B in either graph output and inputs are not correctly entered as a scalar or character. 
2c 95 ‘b’ should be iterated over a larger range to produce the graph
overall: 56


% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% the zero means that when there's no division when there's no population
% and the equation will turn to zero. one is the maximum so the dicition
% rate will stop after it reaches that limit.

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

f=@(a,x) a*x.*(1-x);
interval=(0:0.1:1);
result=f(1,interval); 
figure (1); plot(interval,result);
title('a = 1');

result=f(5,interval);
figure (2); plot(interval,result);
title('a = 5');

result=f(0.3,interval);
figure (3); plot(interval,result);
title('a = 0.3');

% a does have an impact

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 

[timer,xtimes]=ma(1,0.01);

% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

for a=0:0.2:4
for ii=1:200
xo=rand;
x1=a.*xo.*(1-(xo));
for jj=1:500
x2=a.*x1.*(1-(x1));
x3=a.*x2.*(1-(x2));
end
plot(a,x3, '.r');
hold on;
end
end

% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 

% a+b=4
% da/dt=0.25*b*t*v+1
% db/dt=0.25*a*t*v+1

%
% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 
%

solA=ode23(A,[0 4],1);
solB=ode23(B,[0 4],3);
plot(solA.x,solA.y); hold on; plot(solB.x,solB.y);
xlabel('Time'); 
ylabel('Expression');

solA=ode23(A,[0 4],3);
solB=ode23(B,[0 4],1);
plot(solA.x,solA.y); hold on; plot(solB.x,solB.y);
xlabel('Time'); 
ylabel('Expression');

% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 

figure; hold on;
a=1;
for b=0:0.05:3
    polycoeff=[1 -b -a 1];
    rts=roots(polycoeff);
    rts=rts(imag(rts)==0);
    plot(b*ones(length(rts),1),rts,'.r');
end
hold off;
xlabel('b'); 
ylabel('fixed points for v');
