% Define variables and initial guess
syms x y(x) p(x) integrand(y,p)
integrand = sqrt((1+p^2)/(19.62*(1-y)));
first_int = subs(integrand,p,diff(y,x));
second_int = simplify(subs(first_int,y,@(x)-x+1));
initial = vpa(int(second_int,0,1))

% First iteration of functional gradient descent
[new,num] = grad(integrand,-x+1,10);
X = 1*10^(-69):0.005:1;
plot(X,arrayfun(@(x) new(x),X),'r');
title("Gradient Descent vs Optimal");
hold on;
first_descent = num

% Second iteration of functional gradient descent
[new,num] = grad(integrand,new,3);
second_test = num
plot(X,arrayfun(@(x) new(x),X),'g');
hold on;  

% cursed
[new,num] = grad(integrand,new,1.5);
turd_descent = num
plot(X,arrayfun(@(x) new(x),X),'b');
hold on;

% Plot of brachistochrone
syms t x(t) y(t) integrand(y,p)
T = 0:0.005:2.412;
x = @(t) 0.572917*(t-sin(t));
Y = @(t) 0.572917*(cos(t)-1) + 1;
integrand = sqrt((1+p^2)/(19.62*(1-y)));
plot(arrayfun(@(t) x(t),T),arrayfun(@(t) Y(t),T),'k');
hold on;
opt_int = subs(integrand,p,diff(Y,t)/diff(x,t));
second_opt_int = simplify(subs(opt_int,y,Y)*diff(x,t));
optimal = vpaintegral(second_opt_int,0,2.41201)

% Functional gradient descent method
function [eq,value] = grad(integrand,sub,e)
    syms x y(x) p(x) first(p,y) second_int(x) corrector(x)
    % The corrector is a function that is zero valued at the initial/boundary conditions 
    % and positive valued in between. It must also converge to zero faster than the E-L equation 
    % diverges at the initial/boundary conditions (if it diverges at all). This is why the corrector 
    % function has power 3 instead of power 2 in the brachistochrone case.
    corrector = -(x-1)*x^2;
    first = simplify(y - e*corrector*(diff(integrand,y) - diff(diff(integrand,p),x)));
    second = subs(first,p,diff(y,x));
    eq = simplify(subs(second,y,sub));
    first_int = subs(integrand,p,diff(y,x));
    second_int = simplify(subs(first_int,y,eq));
    value = vpaintegral(second_int,0,1);
end

 



