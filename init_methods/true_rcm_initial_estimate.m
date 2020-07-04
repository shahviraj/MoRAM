function x_0 = true_rcm_initial_estimate(A,y_mod,s,R,del,normz,indx)

     [m,n] = size(A);
%     Marg = zeros(1,n); %marginals
%     MShat = zeros(del*s,1); %truncated correlation matrix
%     AShat = zeros(m,del*s); %truncated sensing matrix
%     supp = zeros(1,n); %indicator for initial support Shat
%     y_mod2 = y_mod.^2; %quadratic measurements
% 
%     Marg = ((A'.^2)*(y_mod2))/m; % n x 1
%     [Mg MgS] = sort(Marg,'descend');
%     S0 = MgS(1:del*s); %pick top s-marginals
%     Shat = sort(S0); %store indices in sorted order
%     
%     %supp_recvr(j) = sum(ismember(Shat, z_ind'));
%     AShat = A(:,Shat);
    
    MShat = zeros(n,1); %truncated correlation matrix
    AShat = A;
    
    
    %card_Marg = m;
    %Io = 1:card_Marg;
    
    
    indx_corrected = [indx.idx_0;indx.idx_1;indx.idx_2;indx.idx_3];
    for i = 1:length(indx_corrected)
        ii = indx_corrected(i);
        MShat = MShat + (y_mod(ii))*AShat(ii,:)'; % (s x 1)
    end
    
    MShat = MShat/length(indx_corrected);
    
    %MShat = MShat/(1-(R/2.0)*sqrt(2.0/pi));
%     x_0 = zeros(n,1);
%     x_0(Shat,1)= MShat;
%     x_0 = x_0/(1-(R/2.0)*sqrt(2.0/pi));
    x_0 = MShat;
    x_0 = normz*(x_0/norm(x_0));
end