function eocs = computeEocs (herr)
    % herr = ( (h1,e1),(h2,e2),...,(hm,em) )
    m = length(herr);
    eocs = zeros( 1, m-1 );
    for i = 1:m-1
        eocs(i) = log( herr(i+1,2) / herr(i,2) ) / ...
                  log( herr(i+1,1) / herr(i,1) );
    end
end
