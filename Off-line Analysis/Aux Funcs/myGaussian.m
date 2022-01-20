function f = myGaussian(x,mu,sigma)

p1 = -.5 * ((x - mu)/sigma) .^ 2;
p2 = (sigma * sqrt(2*pi));
f = exp(p1) ./ p2; 


end