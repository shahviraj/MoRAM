function x_s = make_sparse(x,s)
    % make the vector sparse through hard thresholding
    
    n = size(x,1);
    x_abs = abs(x);
    [sorted_x, indx] = sort(x_abs, 'descend');
    
    x_s = zeros(n,1);
    x_s(indx(1:s)) = x(indx(1:s));
    