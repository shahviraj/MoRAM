function p = flip_random(p_in, flip_frac)
    % flip `flip_frac` fraction of values in p_in    
    m = length(p_in);
    k = int64(flip_frac*m);
    p = p_in;
    indx = randperm(m, k);
    p(indx) = 1 - p_in(indx);

end
