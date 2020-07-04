function [x_0,p_refined] = true_rcm_init(A,y_mod,s,pr)
    [m,n] = size(A);
    
    
%     t0 = 0;
%     tl = max(0.0, (pr.R-pr.rho));
%     tu = max(pr.R,pr.rho);
%     tr = pr.R;
    
    % Re-Calculated Measurements
    
    indx.idx_0 = find(y_mod<0.0);
    indx.idx_1 = find(y_mod>=0.0 & y_mod<(pr.R-pr.rho));
    indx.idx_2 = find(y_mod>= pr.rho & y_mod<pr.R);
    indx.idx_3 = find(y_mod>=pr.R);
    %idx_r = m - (indx.idx_0 + indx.idx_1 + indx.idx_2 + indx.idx_3);
    
    y_refined = y_mod;
    p_refined = zeros(m,1);
     
    y_refined(indx.idx_0) = y_mod(indx.idx_0)-pr.R;
    p_refined(indx.idx_0) = 1;
    
    y_refined(indx.idx_1) = y_mod(indx.idx_1);
    p_refined(indx.idx_1) = 0;
    
    y_refined(indx.idx_2) = y_mod(indx.idx_2)-pr.R;
    p_refined(indx.idx_2) = 1;
    
    y_refined(indx.idx_3) = y_mod(indx.idx_3);
    p_refined(indx.idx_3) = 0;
    
    x_0 = true_rcm_initial_estimate(A,y_refined,s,pr.R,pr.del,pr.amp,indx);
    
end

    