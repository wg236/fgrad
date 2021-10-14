% % Brachistochrone Example
% % Define variables and initial guess
syms x y(x) p(x) integrand(y,p)
% integrand = sqrt((1+p^2)/(19.62*(1-y)));
% first_int = subs(integrand,p,diff(y,x));
% second_int = simplify(subs(first_int,y,@(x)-x+1));
% initial = vpa(int(second_int,0,1))
% 
% % First iteration of functional gradient descent
% [new,num] = grad(integrand,-x+1,10);
% X = 1*10^(-69):0.005:1;
% % plot(X,arrayfun(@(x) new(x),X),'r');
% % title("Gradient Descent vs Optimal");
% % hold on;
% first_descent = num
% 
% % Second iteration of functional gradient descent
% [new,num] = grad(integrand,new,3);
% % plot(X,arrayfun(@(x) new(x),X),'g');
% % hold on;  
% second_descent = num
% 
% % cursed
% [new,num] = grad(integrand,new,1.5);
% % plot(X,arrayfun(@(x) new(x),X),'b');
% % hold on;
% turd_descent = num
% 
% % Plot of brachistochrone
% syms t x(t) Y(t) 
% T = 0:0.005:2.41201;
% x = @(t) 0.572917*(t-sin(t));
% Y = @(t) 0.572917*(cos(t)-1) + 1;
% % plot(arrayfun(@(t) x(t),T),arrayfun(@(t) Y(t),T),'k');
% % hold on;
% opt_int = subs(integrand,p,diff(Y,t)/diff(x,t));
% second_opt_int = simplify(subs(opt_int,y,Y(t))*diff(x,t));
% optimal = vpaintegral(second_opt_int,t,0,2.41201)

% Simple Harmonic Oscillator Example
G = 0:5*10^(-3):1;
integrand = p^2 - y^2;
first_int = subs(integrand,p,diff(y,x));
second_int = simplify(subs(first_int,y,@(x)x));
initial = vpa(int(second_int,0,1))
[new,num] = grad(integrand,x,0.3);
for i=1:27
    [new,num] = grad(integrand,new,0.01);
end
% plot(G,arrayfun(@(x) sin(x)/sin(1), G),'k');
% hold on;
% plot(G,arrayfun(@(x) new(x),G),'r');
% hold on;
el = diff(integrand,y) - diff(diff(integrand,p),x);
el = subs(el,p,diff(y,x));
el = subs(el,y,new);
average_dev_1 = vpaintegral(abs(el),x,0,1)
plot(G,arrayfun(@(x) el(x),G),'r');
hold on;
[new,num] = grad(integrand,x,0.3);
for i=1:17
    [new,num] = grad(integrand,new,0.02);
end
el = diff(integrand,y) - diff(diff(integrand,p),x);
el = subs(el,p,diff(y,x));
el = subs(el,y,new);
average_dev_2 = vpaintegral(abs(el),x,0,1)
plot(G,arrayfun(@(x) el(x),G),'g');
ylim([-0.6 0.6]);
yticks(-0.6:0.2:0.6);

opt_int = subs(integrand,p,cos(x)/sin(1));
second_opt_int = simplify(subs(opt_int,y,sin(x)/sin(1)));
optimal = vpaintegral(second_opt_int,x,0,1)

% Functional gradient descent method
function [eq,value] = grad(integrand,sub,e)
    syms x y(x) p(x) first(p,y) second_int(x) corrector(x)
    % The corrector is a function that is zero valued at the initial/boundary conditions 
    % and positive valued in between. It must also converge to zero faster than the E-L equation 
    % diverges at the initial/boundary conditions (if it diverges at all) 
    corrector = -(x-1)*x;
    first = simplify(y-e*corrector*(diff(integrand,y) - diff(diff(integrand,p),x)));
    second = subs(first,p,diff(y,x));
    eq = simplify(subs(second,y,sub));
    first_int = subs(integrand,p,diff(y,x));
    second_int = simplify(subs(first_int,y,eq));
    value = vpaintegral(second_int,0,1);
end

 



