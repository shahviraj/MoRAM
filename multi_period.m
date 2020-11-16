% This file is for experimenting with various new ideas : mainly for EURASI
% P response
clear all;
clc;
tic;
rng ('shuffle')
pr = struct;
%Fixed parameters
pr.n = 1000; %length of the input signal
pr.b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
pr.tol1 = 1e-8; %error tolerance for measurements
pr.tol2 = 1e-8;
pr.max_iter = 15;
pr.flip_frac = 0.01;
pr.R = 3; %period of the modulo function
pr.rho = 3;%spread of the true measurements, y =A*z
pr.del = 1; %truncation factor for supp estimation
pr.spgl_opts = spgSetParms('verbosity',0);
%Tuned parameters
%pr.mspan1 = [100:100:2000];
%pr.mspan2 = [600:100:1000];
%pr.mspan=[pr.mspan1,pr.mspan2];
pr.m =500;
pr.num_trials = 1;
pr.s = 18; % sparsity
pr.amp = 1; %amplification factor 
pr.del_p = 0.15; % ps = del*m (sparsity pertaining to error in p)

m = pr.m;

s = pr.s;
ps = floor(pr.del_p*m);
%Generate the ground truth signal
[z,z_ind] =  generate_signal(pr.n,s,pr.b, pr.amp);

%Generate the measurements: y=mod(Ax,R)
[y_mod,y_p, A] = multi_period_modulo_measure_signal(m,z,pr.R);

% Calculate the initial estimate
p_refined = simple_rcm_init(A,y_mod,s,pr);

% we give initial p as calculated from simple rcm init
p = p_refined;
fprintf('\n#iter\t\t|y-Ax|\t\t|x-z|\trecovery_prob\n')

for t=1:pr.max_iter
    
    % removed x from the input of the function below
    if t < 16
        [x,delta_p] = mod_spgl1_bp(y_mod,p,A, pr.R, pr.spgl_opts); % SPGL1 implementation -- faster
    end
%     if t > 3
%         [x,delta_p] = mod_cosamp(y_mod,p,A,x,pr.R,s,ps);
%     end
    
    p = (-sign(A*x)+1)/2;

    %p = p - delta_p;
%     if t > 2
%         p = flip_random(p,pr.flip_frac);
%     end
    err_hist(t+1,1) = norm((y_mod-y_p*pr.R)-(A*x))/norm(y_mod-y_p*pr.R);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    recovery_prob = nnz(~(p-y_p));
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\t\t%2.0f\n',t,err_hist(t+1,1),err_hist(t+1,2),recovery_prob)
%     if (err_hist(t+1,1) < pr.tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<pr.tol2)
%         break;
%     end
end
% Relative reconstruction error
reconst_err = norm(x-z)/norm(z)