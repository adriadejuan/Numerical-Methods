function y = evolve(phi,f,Df, t0,y0, T,N)
    h = T/N;
    y = zeros([ N+1 length(y0) ]);
    y(1) = y0;
    t = 0;
    for i=1:N
        y(i+1) = phi(f,Df, t,y(i), h);
        t = t+h;
    end
end
