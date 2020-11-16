clear all;
clc;
tic;
%rng ('shuffle')
pr = struct;
%Fixed parameters
pr.n = 1000; %length of the input signal
pr.b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
pr.tol1 = 1e-5; %error tolerance for measurements
pr.tol2 = 1e-7;
pr.max_iter = 5;
pr.r_span = 3:0.1:4.0; %period of the modulo function
pr.rho = 3;%spread of the true measurements, y =A*z
pr.del = 1; %truncation factor for supp estimation
pr.spgl_opts = spgSetParms('verbosity',0);
%Tuned parameters
%pr.mspan1 = [100:100:2000];
%pr.mspan2 = [600:100:1000];
%pr.mspan=[pr.mspan1,pr.mspan2];
pr.m=400;
pr.num_trials = 20;
pr.s_span = 12:12:12; % sparsity
pr.amp = 1; %amplification factor 
pr.del_p = 0.005; % ps = del*m (sparsity pertaining to error in p)
pr.method = 'justice-pursuit';
pr.init_method = 'simple-rcm';
pr.svd_opt = 'svd';
pr.plot_method = 'mean-error';

all_rho = zeros(length(pr.mspan),length(pr.s_span),pr.num_trials);

err = zeros(length(pr.mspan),length(pr.s_span),pr.num_trials);
supp_recvr=zeros(length(pr.mspan),length(pr.s_span));
init_err=zeros(length(pr.mspan),length(pr.s_span),pr.num_trials);
reconst_err=zeros(length(pr.mspan),length(pr.s_span));

for j = 1:length(pr.r_span)
    m = pr.m;
    R = pr.r_span(j);
    for k = 1:length(pr.s_span)
        s = pr.s_span(k);
        
        for l = 1:pr.num_trials
            %Generate the ground truth signal
            [z,z_ind] =  generate_signal(pr.n,s,pr.b, pr.amp);

            %Generate the measurements: y=mod(Ax,R)
            [y_mod, y, y_p, A] = modulo_measure_signal(m,z,pr.R);
            all_rho(j,k,l) = max(abs(y));
            switch pr.init_method
                case 'copram'
                    x_0 = copram_init(y_mod,A,s);
                case 'moram'
                    x_0 = moram_init(A,y_mod,s,pr.R,pr.del,pr.amp);
    %             case 'raf'
    %                 x_0 = raf_init();
%                 case 'rcm' %Re-Calculated Measurements
%                     [x_0,p_refined,idx] = rcm_init(A,y_mod,s,pr);
                case 'true-rcm'
                    [x_0,p_refined] = true_rcm_init(A,y_mod,s,pr);
                    %extra line added to make x_0 sparse
                    x_0 = make_sparse(x_0,s);
                case 'simple-rcm'
                    p_refined = simple_rcm_init(A,y_mod,pr);
                    %extra line added to make x_0 sparse
                    %x_0 = make_sparse(x_0,s);
            end


            fprintf('\n#iter\t\t|y-Ax|\t\t|x-z|\trecovery_prob\n')
            
            p = p_refined;
            for t=1:pr.max_iter
                switch pr.method
                    case 'cosamp'
                        x = cosamp((y_mod-pr.R*p)/sqrt(m), A/sqrt(m),s,100,x); %Its = 100

                        y_eff = (y_mod-pr.R*p)/sqrt(m);
                        A_eff = A/sqrt(m);
                        disp('error in CoSaMP output')
                        norm(y_eff-(A_eff*x))/norm(y_eff)

                    case 'justice-pursuit'
                        %[x,delta_p] = mod_l1_bp(y_mod,p,A,x,pr.R); % l1 -magic implementation -- slow
                        
                        [x,delta_p] = mod_spgl1_bp(y_mod,p,A,pr.R, pr.spgl_opts); % SPGL1 implementation -- faster
                        
                    case 'basis-pursuit'
                        x = l1eq_pd(x,A/sqrt(m), [], (y_mod-pr.R*p)/sqrt(m)); % l1-magic implementation -- slow
                        
                    case 'robust-cosamp'
                        [x,delta_p] = mod_cosamp(y_mod,p,A,x,pr.R,s,ps);
                end
                p = (-sign(A*x)+1)/2;
           
                err_hist(t+1,1) = norm((y_mod-y_p*pr.R)-(A*x))/norm(y_mod-y_p*pr.R);
                err_hist(t+1,2) = norm(x-z)/norm(z);
                recovery_prob = nnz(~(p-y_p));
                fprintf('\n%d\t\t%2.8f\t\t%2.4f\t\t%2.0f\n',t,err_hist(t+1,1),err_hist(t+1,2),recovery_prob)
                if (err_hist(t+1,1) < pr.tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<pr.tol2)
                    break;
                end
            end
            % Relative reconstruction error
            reconst_err(j,k,l) = norm(x-z)/norm(z);
        end
        
    end
    
end
toc
% p_err = y_p - p;
% p_err_idx = find(p_err~=0);
% y_true = A*z;
% disp('value of y corresponding to the errors are:')
% disp(y_true(p_err_idx))
% 
% 
% y_sorted = sort(y_true, 'ComparisonMethod','abs');
% y_sorted(1:length(p_err_idx))

if ~exist('./results', 'dir')
       mkdir('./results')
end

meanerr = squeeze(mean(reconst_err,3));
for i = 1:size(meanerr,2)
    fprintf('\nSparsity: %d\t\t\n', pr.s_span(i));
    disp("Measurements");
    disp(pr.mspan); 
    disp("Mean errors");
    disp(meanerr(1:length(pr.mspan),i));
end
% construct_subplots(reconst_err,pr, ['rconst_',pr.init_method,'_amp_',num2str(pr.amp),'_r_',num2str(pr.R),'_s_',...
%     num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
%     '_',num2str(pr.mspan(end)),'_',pr.method,'_num_trials_',num2str(pr.num_trials)],pr.plot_method,1);

% construct_plot(init_err,pr,['init_','r_',num2str(pr.R),'_s_',...
%     num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
%     '_',num2str(pr.mspan(end)),'_',pr.method],pr.method);