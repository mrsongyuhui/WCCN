clc;
clear all;

for t = 1:50
    % Parameter
    cell_size = 500;
    M = 2;
    N = 2.^(2:2:20);
    K = 38;
    tau = 36;
    P = 4;
    S = 36;
    SNR_db = 111;
    SNR_linear = 10^(SNR_db/10);
    % Simulation
    for i = 1:size(N,2)
        [rate_CF_C_OP(i,t,:),rate_CF_C_OR(i,t,:),rate_CF_D_OP(i,t,:),rate_CF_D_OR(i,t,:),rate_perfect(i,t,:)] = RATE(cell_size,M,N(i),K,tau,P,S,SNR_linear,'CF','OFDM');
    end
end

figure(1);
plot(log2(N),mean(rate_perfect,2),'LineWidth',2,'Color','#77AC30','LineStyle','-','Marker','o'); hold on;
plot(log2(N),mean(rate_CF_C_OP,2),'LineWidth',2,'Color','#0072BD','LineStyle','-.'); hold on;
plot(log2(N),mean(rate_CF_C_OR,2),'LineWidth',2,'Color','#4DBEEE','LineStyle',':');
plot(log2(N),mean(rate_CF_D_OP,2),'LineWidth',2,'Color','#A2142F','LineStyle','-'); hold on;
plot(log2(N),mean(rate_CF_D_OR,2),'LineWidth',2,'Color','#D95319','LineStyle','--');
legend('perfect','co\_op','co\_fea','non\_op','non\_fea','interpreter','latex','FontSize', 20, 'FontName', 'Times New Roman')
xlabel('$\log_2\left(N\right)$','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Achievable rate (bits/s/Hz)','interpreter','latex', 'FontSize', 20, 'FontName', 'Times New Roman');

function [rate_C_OP,rate_C_OR,rate_D_OP,rate_D_OR,rate_perfect] = RATE(cell_size,M,N,K,tau,P,S,SNR_linear,cell_mode,pilot_transmission)
%% Generate channel
h = cellfun(@(~)cellfun(@(~)(randn(N,P)+1i*randn(N,P))/sqrt(2),cell(K,1),'UniformOutput',false),cell(M,1),'UniformOutput',false);
switch cell_mode
    case 'CF'
        coordinate = Hexagon_Random(cell_size, 50, 7, K/2, 2);
        base_location = cell2mat(cellfun(@(a)cell2mat(a(2,1)),coordinate,'UniformOutput',false)).';
        user_location = cell2mat(cellfun(@(a)cell2mat(a(3,1)),coordinate,'UniformOutput',false)).';
        distance = abs(base_location-user_location.');
        beta = 10.^(-(128.1+37.6*log10((distance/1000)))/10);
    case 'CO'
        beta = (1-0.6).*rand(M,K) + 0.6;
end

%% Transmission
% Transmit pilot
zc_lengtth = 2*floor(tau/2)+1;
potential_primes = primes(zc_lengtth);
final_primes = [1 potential_primes(gcd(zc_lengtth,potential_primes)==1)];
switch pilot_transmission
    case 'OFDM'
        x_p_time_domain_all = cellfun(@(b)b(:,1:P).',cellfun(@(a)pinv(dftmtx(tau))*diag(dftmtx(tau)*a)*dftmtx(tau),mat2cell(cell2mat(cellfun(@(d)d(1:tau,:),cellfun(@(c)pinv(dftmtx(zc_lengtth))*diag(dftmtx(zc_lengtth)*c)*dftmtx(zc_lengtth),mat2cell(exp(-1i*pi*linspace(0,zc_lengtth-1,zc_lengtth)'.*(linspace(0,zc_lengtth-1,zc_lengtth)'+1)*final_primes/zc_lengtth),zc_lengtth,ones(1,size(final_primes,2))),'UniformOutput',false),'UniformOutput',false)),tau,ones(1,(size(final_primes,2))*zc_lengtth)),'UniformOutput',false),'UniformOutput',false)';
    case 'Non_OFDM'
        x_p_time_domain_all = cellfun(@(a)toeplitz(a,[a(1) zeros(1,P-1)]).',mat2cell(cell2mat(cellfun(@(d)d(1:tau,:),cellfun(@(c)pinv(dftmtx(zc_lengtth))*diag(dftmtx(zc_lengtth)*c)*dftmtx(zc_lengtth),mat2cell(exp(-1i*pi*linspace(0,zc_lengtth-1,zc_lengtth)'.*(linspace(0,zc_lengtth-1,zc_lengtth)'+1)*final_primes/zc_lengtth),zc_lengtth,ones(1,size(final_primes,2))),'UniformOutput',false),'UniformOutput',false)),tau,ones(1,(size(final_primes,2))*zc_lengtth)),'UniformOutput',false)';
end
x_p_time_domain = x_p_time_domain_all(1:K);
y_p_clean=cellfun(@(d,e)sqrt(SNR_linear)*sum(cat(3,e{:}).*reshape(d,1,1,[]),3),mat2cell(sqrt(beta),ones(M,1),K),cellfun(@(c)cellfun(@(a,b)a*b,c,x_p_time_domain,'UniformOutput',false),h,'UniformOutput',false),'UniformOutput',false);
y_p = cellfun(@(a)a+(randn(N,tau)+1i*randn(N,tau))/sqrt(2),y_p_clean,'UniformOutput',false);
% DFT matrix
F_all = dftmtx(S);
F = F_all(1:P,:);
data = sign(randn(1,S)); all_data = repmat(data',1,K);%.*[ones(S,1),zeros(S,K-1)];
% Transmit data
y_d_clean=cellfun(@(d,e)sqrt(SNR_linear)*sum(cat(3,e{:}).*reshape(d,1,1,[]),3),mat2cell(sqrt(beta),ones(M,1),K),cellfun(@(c)cellfun(@(a)a*F.*repmat(data,N,1),c,'UniformOutput',false),h,'UniformOutput',false),'UniformOutput',false);
y_d = cellfun(@(a)a+(randn(N,S)+1i*randn(N,S))/sqrt(2),y_d_clean,'UniformOutput',false);

%% Estimator and Detection
Q = cell2mat(cellfun(@(e)cat(1,e{:}),{cellfun(@(d)cat(3,d{:}),cellfun(@(a)cellfun(@(b,c)b.'*conj(F)*N*c*sqrt(SNR_linear),x_p_time_domain',num2cell(a),'UniformOutput',false),mat2cell(beta,ones(M,1),K),'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false));
% 'estimate_centralized_original'
% Estimator
estimator_all = cell2mat(cellfun(@(b)cat(3,b{:}),{cellfun(@(a)calculate(Q(:,a,:)),num2cell(linspace(1,S,S)),'UniformOutput',false)},'UniformOutput',false));
estimator_seperatecentralized_original = cellfun(@(a)cellfun(@(b)squeeze(b),mat2cell(a,tau,ones(K,1),S),'UniformOutput',false),mat2cell(estimator_all,ones(M,1)*tau,K,S),'UniformOutput',false);
estimated_channelcentralized_original = cellfun(@(a,b)cellfun(@(c)a*c,b,'UniformOutput',false),y_p,estimator_seperatecentralized_original,'UniformOutput',false);
% Detection
detected_data = sum(cell2mat(cellfun(@(e)cat(3,e{:}),{cellfun(@(d)cell2mat(d),cellfun(@(a,b)cellfun(@(c)diag(c'*a),b,'UniformOutput',false),y_d,estimated_channelcentralized_original,'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
scale = sum(cell2mat(cellfun(@(g)cat(3,g{:}),{cellfun(@(f)cell2mat(f),cellfun(@(a,b)cellfun(@(c,d,e)diag(c'*e'*F*N*d*SNR_linear),a,num2cell(b),x_p_time_domain','UniformOutput',false),estimator_seperatecentralized_original,mat2cell(beta,ones(M,1),K),'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
rate_C_OR = mean(log2(1+ones(S,K)./abs(detected_data./scale-all_data).^2),'all');

% 'estimate_distributed_original'
% Estimator
estimator_seperate_distributed_original = cellfun(@(d)mat2cell(d,tau,ones(1,K),S),cellfun(@(c)cat(3,c{:}),cellfun(@(b)cellfun(@(a)calculate(b(:,a,:)),num2cell(linspace(1,S,S)),'UniformOutput',false),mat2cell(Q,ones(M,1)*tau,S,K),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
estimated_channel_distributed_original = cellfun(@(a,b)cellfun(@(c)a*squeeze(c),b,'UniformOutput',false),y_p,estimator_seperate_distributed_original,'UniformOutput',false);
% Detection
detected_data = sum(cell2mat(cellfun(@(e)cat(3,e{:}),{cellfun(@(d)cell2mat(d),cellfun(@(a,b)cellfun(@(c)diag(c'*a),b,'UniformOutput',false),y_d,estimated_channel_distributed_original,'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
scale = sum(cell2mat(cellfun(@(g)cat(3,g{:}),{cellfun(@(f)cell2mat(f),cellfun(@(a,b)cellfun(@(c,d,e)diag(squeeze(c)'*e'*F*N*d*SNR_linear),a,num2cell(b),x_p_time_domain','UniformOutput',false),estimator_seperate_distributed_original,mat2cell(beta,ones(M,1),K),'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
rate_D_OR = mean(log2(1+ones(S,K)./abs(detected_data./scale-all_data).^2),'all');

% 'estimate_centralized_optimized'
estimator_KSM = cellfun(@(k)cellfun(@(s)mat2cell(cell2mat(cellfun(@(d)...
    rayleigh_solver(d*d'*SNR_linear,...
    cell2mat(cellfun(@(m,n)numerator(m,n,k,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F(:,s)),num2cell(repmat(linspace(1,M,M)',1,M)),num2cell(repmat(linspace(1,M,M),M,1)),'UniformOutput',false))...
    -d*d'*SNR_linear...
    +cell2mat(cellfun(@(m,n)denominator(m,n,k,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F(:,s)),num2cell(repmat(linspace(1,M,M)',1,M)),num2cell(repmat(linspace(1,M,M),M,1)),'UniformOutput',false))),...
    {cell2mat(cellfun(@(m)(x_p_time_domain{k}'*beta(m,k)*N*F(:,s)),num2cell(linspace(1,M,M)),'UniformOutput',false)')}...
    ,'UniformOutput',false)),ones(M,1)*tau,1)',num2cell(linspace(1,S,S)),'UniformOutput',false),num2cell(linspace(1,K,K)),'UniformOutput',false);

estimator_matrix = cell2mat(cellfun(@(d)cat(4,d{:}),{cellfun(@(c)cell2mat(cellfun(@(b)cat(3,b{:}),{cellfun(@(a)cell2mat(a),c,'UniformOutput',false)},'UniformOutput',false)),estimator_KSM,'UniformOutput',false)},'UniformOutput',false));
estimator_seperate_distributed_optimized = cellfun(@(a)squeeze(mat2cell(squeeze(a),tau,S,ones(1,K))).',mat2cell(estimator_matrix,tau,ones(M,1),S,K),'UniformOutput',false).';
estimated_channel_distributed_optimized = cellfun(@(a,b)cellfun(@(c)a*c,b,'UniformOutput',false),y_p,estimator_seperate_distributed_optimized,'UniformOutput',false);
% Detection
detected_data = sum(cell2mat(cellfun(@(e)cat(3,e{:}),{cellfun(@(d)cell2mat(d),cellfun(@(a,b)cellfun(@(c)diag(c'*a),b,'UniformOutput',false),y_d,estimated_channel_distributed_optimized,'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
scale = sum(cell2mat(cellfun(@(g)cat(3,g{:}),{cellfun(@(f)cell2mat(f),cellfun(@(a,b)cellfun(@(c,d,e)diag(squeeze(c)'*e'*F*N*d*SNR_linear),a,num2cell(b),x_p_time_domain','UniformOutput',false),estimator_seperate_distributed_optimized,mat2cell(beta,ones(M,1),K),'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
rate_C_OP = mean(log2(1+ones(S,K)./abs(detected_data./scale-all_data).^2),'all');

% 'estimate_distributed_optimized'
% Estimator
estimator_KSM = cellfun(@(k)cellfun(@(s)cellfun(@(m)rayleigh_solver(numerator(m,m,k,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F(:,s)),denominator(m,m,k,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F(:,s))),num2cell(linspace(1,M,M)),'UniformOutput',false),num2cell(linspace(1,S,S)),'UniformOutput',false),num2cell(linspace(1,K,K)),'UniformOutput',false);
estimator_matrix = cell2mat(cellfun(@(d)cat(4,d{:}),{cellfun(@(c)cell2mat(cellfun(@(b)cat(3,b{:}),{cellfun(@(a)cell2mat(a),c,'UniformOutput',false)},'UniformOutput',false)),estimator_KSM,'UniformOutput',false)},'UniformOutput',false));
estimator_seperate_distributed_optimized = cellfun(@(a)squeeze(mat2cell(squeeze(a),tau,S,ones(1,K))).',mat2cell(estimator_matrix,tau,ones(M,1),S,K),'UniformOutput',false).';
estimated_channel_distributed_optimized = cellfun(@(a,b)cellfun(@(c)a*c,b,'UniformOutput',false),y_p,estimator_seperate_distributed_optimized,'UniformOutput',false);
% Detection
detected_data = sum(cell2mat(cellfun(@(e)cat(3,e{:}),{cellfun(@(d)cell2mat(d),cellfun(@(a,b)cellfun(@(c)diag(c'*a),b,'UniformOutput',false),y_d,estimated_channel_distributed_optimized,'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
scale = sum(cell2mat(cellfun(@(g)cat(3,g{:}),{cellfun(@(f)cell2mat(f),cellfun(@(a,b)cellfun(@(c,d,e)diag(squeeze(c)'*e'*F*N*d*SNR_linear),a,num2cell(b),x_p_time_domain','UniformOutput',false),estimator_seperate_distributed_optimized,mat2cell(beta,ones(M,1),K),'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
rate_D_OP = mean(log2(1+ones(S,K)./abs(detected_data./scale-all_data).^2),'all');
% 'perfect'
g = cellfun(@(a,b)cellfun(@(c,d)c*d*F,num2cell(a'),b,'UniformOutput',false)',mat2cell(sqrt(beta),ones(M,1),K),h,'UniformOutput',false);
detected_data = sum(cell2mat(cellfun(@(e)cat(3,e{:}),{cellfun(@(d)cell2mat(d),cellfun(@(a,b)cellfun(@(c)diag(c'*a),b,'UniformOutput',false),y_d,g,'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
scale = sqrt(SNR_linear)*sum(cell2mat(cellfun(@(e)cat(3,e{:}),{cellfun(@(d)cell2mat(d),cellfun(@(a)cellfun(@(b)diag(b'*b),a,'UniformOutput',false),g,'UniformOutput',false),'UniformOutput',false)},'UniformOutput',false)),3);
rate_perfect = mean(log2(1+ones(S,K)./abs(detected_data./scale-all_data).^2),'all');
end

function [result] = rayleigh_solver(A,B)
[V,D] = eig(A,B,'vector');
[~,I] = max(D);
result = V(:,I);
end


function [result] = eigenvalue(A,B)
[~,D] = eig(A,B,'vector');
result = max(D);
end

function [result] = numerator(m,n,k,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F)
switch m==n
    case 1
        result = ((x_p_time_domain{k}'*beta(m,k)*SNR_linear*N*F)*(F'*beta(m,k)*N*x_p_time_domain{k}))+...
            cell2mat(cellfun(@(a)sum(cat(3,a{:}),3),{cellfun(@(j)x_p_time_domain{j}'*beta(m,k)*SNR_linear*P*beta(m,j)*N*x_p_time_domain{j},num2cell(setdiff(linspace(1,K,K),k)),'UniformOutput',false)},'UniformOutput',false))+...
            beta(m,k)*P*N*eye(tau);
    case 0
        result = ((x_p_time_domain{k}'*beta(m,k)*SNR_linear*N*F)*(F'*beta(n,k)*N*x_p_time_domain{k}));
end
end

function [result] = denominator(m,n,k,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F)
switch m==n
    case {0,1}
        result = cell2mat(cellfun(@(b)sum(cat(3,b{:}),3),{cellfun(@(a)numerator(m,n,a,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F),num2cell(setdiff(linspace(1,K,K),k)),'UniformOutput',false)},'UniformOutput',false))+...
            cell2mat(cellfun(@(b)sum(cat(3,b{:}),3),{cellfun(@(a)x_p_time_domain{a}'*beta(m,a)*N*x_p_time_domain{a},num2cell(linspace(1,K,K)),'UniformOutput',false)},'UniformOutput',false))+...
            N*eye(tau)/SNR_linear;
    case 2
        result = cell2mat(cellfun(@(b)sum(cat(3,b{:}),3),{cellfun(@(a)numerator(m,n,a,N,K,tau,P,SNR_linear,x_p_time_domain,beta,F),num2cell(setdiff(linspace(1,K,K),k)),'UniformOutput',false)},'UniformOutput',false));
end
end

function [result] = calculate(A)
[U,D,V] = svd(squeeze(A).');
result = pinv(D*V')*U';
end
