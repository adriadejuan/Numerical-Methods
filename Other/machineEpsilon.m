% COMPUTE MACHINE EPSILON %
% Machine Epsilon is the lowest real number u of the form
% u = 2^-i for i > 0 which satisfies 1 + u > 1

eps = 1;

while 1.0 + eps > 1.0
    eps = eps /2;
end
