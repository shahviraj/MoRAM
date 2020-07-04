function x_0 = raf_init(A,y_mod,s,pr)

    [m,n] = size(A);
    Marg = zeros(1,n); %marginals
    MShat_neg = zeros(pr.del*s); %truncated correlation matrix
    MShat_pos = zeros(pr.del*s); %truncated correlation matrix
    AShat = zeros(m,pr.del*s); %truncated sensing matrix
    supp = zeros(1,n); %indicator for initial support Shat
    y_mod2 = y_mod.^2; %quadratic measurements
    neg_c = 0;
    pos_c = 0;
    
    
    Marg = ((A'.^2)*(y_mod2))/m; % n x 1
    [Mg MgS] = sort(Marg,'descend');
    S0 = MgS(1:pr.del*s); %pick top s-marginals
    Shat = sort(S0); %store indices in sorted order
    
    %supp_recvr(j) = sum(ismember(Shat, z_ind'));
    AShat = A(:,Shat);

    
    %calculate r_hat
    r_hat = sqrt(mean(y_mod2));
    
    
    card_Marg = m;
    Io = 1:card_Marg;
  
    for i = 1:card_Marg
        ii = Io(i);
        if y_mod2(ii)<(r_hat^2)/2
            MShat_neg = MShat_neg + (pr.l_neg)*(AShat(ii,:)'*AShat(ii,:));
            neg_c = neg_c+1;
        else
            MShat_pos = MShat_pos + (pr.l_pos)*(AShat(ii,:)'*AShat(ii,:));
            pos_c = pos_c+1;
        end
    end
    
    MShat = MShat_neg/neg_c + MShat_pos/pos_c; %(s x s)

%     svd_opt = 'svd'; %more accurate, but slower for larger dimensions
%     svd_opt = 'power'; %approximate, faster for larger dimensions

    switch pr.svd_opt
        case 'svd'
            [u,sigma,v] = svd(MShat);
            v1 = u(:,1); %top singular vector of MShat, normalized - s x 1
        case 'power'
            v1 = svd_power(MShat);
    end

    v = zeros(n,1);
    v(Shat,1) = v1;
    x_0 = r_hat*v; %ensures that the energy/norm of the initial estimate is close to actual
    %x_0 = pr.amp*v;
%     x_0 = zeros(n,1);
%     x_0(Shat,1)= MShat;
%     x_0 = x_0/(1-(R/2.0)*sqrt(2.0/pi));
%     
%     x_0 = normz*(x_0/norm(x_0));
end

    