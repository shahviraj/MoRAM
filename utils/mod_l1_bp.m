function [x,delta_p] = mod_l1_bp(y_mod,p, A, x, R)

n= length(x);
m = length(y_mod);

%Aps = R*eye(m);
Aps = eye(m);
aug_A = [A/sqrt(m), Aps];

delta_p = zeros(m,1);
aug_x = [x;delta_p];

aug_x = l1eq_pd(aug_x,aug_A, [], (y_mod-R*p)/sqrt(m));
%aug_x = cosamp((y_mod-R*p), aug_A, s+ps ,100, aug_x);
x= aug_x(1:n);
delta_p = aug_x(n+1:end);



