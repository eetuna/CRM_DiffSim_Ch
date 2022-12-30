function Df = jacobian_forward_difference(f,x,h)
% Find a forward difference approximation the Jacobian of a function 
% f:R^n->R^m.
% f and x are column vectors.
% h is the step size.

f_0 = f(x);
n = length(x);
Df = zeros(length(f_0),n);
X = x*ones(1,n) + h*eye(n);
for i = 1:n
    Df(:,i) = f(X(:,i)) - f_0;
end
Df = Df/h;

end
