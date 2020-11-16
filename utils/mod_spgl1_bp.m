function [x,delta_p] = mod_spgl1_bp(y_mod,p, A, R,opts)
% x is removed from the input itself
n = size(A,2);
m = length(y_mod);

%Aps = R*eye(m);
Aps = eye(m);
aug_A = [A/sqrt(m), Aps];

aug_x = spg_bp(aug_A,(y_mod-R*p)/sqrt(m),opts);
%aug_x = cosamp((y_mod-R*p), aug_A, s+ps ,100, aug_x);
x= aug_x(1:n);
delta_p = aug_x(n+1:end);



