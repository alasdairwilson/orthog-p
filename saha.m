%setup constants

% ga = 2;
% gi = 2;
% ge = 1;
 e =  1.602e-19;
 vi =  13.6;
% n = 1e20;
k = 1.38e-23;
m = 9.11e-31;
h = 6.63e-34;
n = 1e20;
pi=3.14159;

%set coefficients in saha eqn for simplicity
a = (1/n)*(2*pi*m*k/(h^2))^1.5;
b = -vi*e/k;


%The following equations can be used to calculate either the ionisation
%fraction as a function of temperature or the temperature as a function of
%the ionisation fraction.

%set T
T = [9004.6]
%calculates LHS of saha
y = a*T^1.5*exp((b)/T);
%calculates ionisation fraction
frac = 1/2*(sqrt(y)*sqrt(y+4)-y)

%set ion frac
frac = 0.5

%calculates LHS of saha
y =  (frac^2/(1-frac));


%calculates temperature from fraction.
T = -(2*b)/(3*lambertw(-(2*a^(2/3)*b)/(3*y^(2/3))))
