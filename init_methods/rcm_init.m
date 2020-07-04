function [x_0,p_refined,idx_r] = rcm_init(A,y_mod,s,pr)
    [m,n] = size(A);
    
    t0 = pr.R;
    t1 = 0.75*pr.R;
    t2 = 0.25*pr.R;
    
    % Re-Calculated Measurements
    
    indx.idx_0 = find(y_mod>=t0);
    indx.idx_1 = find(y_mod>=t1 & y_mod<t0);
    indx.idx_2 = find(y_mod>=t2 & y_mod<t1);
    indx.idx_3 = find(y_mod<t2);
    idx_r = indx.idx_2;
    
    y_refined = zeros(m,1);
    p_refined = zeros(m,1);
     
    y_refined(indx.idx_0) = y_mod(indx.idx_0);
    p_refined(indx.idx_0) = 0;
    
    y_refined(indx.idx_1) = y_mod(indx.idx_1)-pr.R;
    p_refined(indx.idx_1) = 1;
    
    y_refined(indx.idx_2) = y_mod(indx.idx_2);
    p_refined(indx.idx_2) = 0;
    
    y_refined(indx.idx_3) = y_mod(indx.idx_3);
    p_refined(indx.idx_3) = 0;
    
    x_0 = rcm_initial_estimate(A,y_refined,s,pr.R,pr.del,pr.amp,indx);
    
end

    