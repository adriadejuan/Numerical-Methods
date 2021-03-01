function [x,k] = newton(F,DF,x0,eps,K)
    k  = 0;
    x  = x0;
    Fx = F(x);
    while dot(Fx,Fx) > eps*eps && k<K
        x  = x - (DF(x)\Fx')';
        Fx = F(x);
        k  = k+1;
    end
end
