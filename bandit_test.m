function cumulR_soft = bandit_test(n_ban,T,alpha,beta)
%% Agent parameters
% bandit(4,1000,0.1,2)
close all;
rng(3)
%% bandit parameters
% key is expected reward E(i) = P(i) x R(i)
% P = [1 1 1 1];              % deterministic bandit
% R = [0.4 0.7 0.1 0.5];      % variable rewards

P = ones(1,n_ban);
% R = int32(randi([0,1],[1,n_ban])) Question: 0 or 1??
% R = rand(1,n_ban)
R = [0.4 0.7 0.1 0.5];
labels = num2str(num2str([R].','Reward = %.2f'));
% P = [0.4 0.7 0.1 0.5];      % stochastic bandits
% R = [1 1 1 1];              % identical rewards (== same expected reward as above)

N = numel(P);

%% learn bandits

R_soft = zeros(1,T);
A_soft = zeros(1,T);

% initialise learnt values
Q_soft = zeros(T,N);
Q_soft(1,:) = rand(1,N) * 0.1;  % initialise learnt values to some small random number
% %%constant beta

% for t = 1:T
%     % choose action (pull arm of bandit) based on current estimated value of each bandit
%     P_action_soft = exp(Q_soft(t,:)*beta) ./ sum(exp(Q_soft(t,:)*beta)); % softmax P(action)
%     A_soft(t) = find(rand <= cumsum(P_action_soft),1);  % choose action based on those probabilities
%     
%     % was action rewarded?
%     R_soft(t) = (rand <= P(A_soft(t))) * R(A_soft(t));  
%     % Update values of actions giiven reward:
%     % (1) all actions not chosen keep same value estimates
%     Q_soft(t+1,:) = Q_soft(t,:);
%     
%     % (2) update value estimate of chosen action: Q(new) = Q(old) + alpha * (R - Q(old))
%     Q_soft(t+1,A_soft(t)) = Q_soft(t,A_soft(t)) + alpha * (R_soft(t) - Q_soft(t,A_soft(t)));
%     
% end
%% dynamic beta
% dynamic_exp = exp_beta(2,100); %exp_beta(A,tau)
dynamic_lin = lin_beta(0.001);

for t = 1:T
    % choose action (pull arm of bandit) based on current estimated value of each bandit
%     P_action_soft = exp(Q_soft(t,:)*dynamic_exp(t)) ./ sum(exp(Q_soft(t,:)*dynamic_exp(t))); % softmax P(action)
    P_action_soft = exp(Q_soft(t,:)*dynamic_lin(t)) ./ sum(exp(Q_soft(t,:)*dynamic_lin(t))); % softmax P(action)
    A_soft(t) = find(rand <= cumsum(P_action_soft),1);  % choose action based on those probabilities
    
    % was action rewarded?
    R_soft(t) = (rand <= P(A_soft(t))) * R(A_soft(t));  
    % Update values of actions giiven reward:
    % (1) all actions not chosen keep same value estimates
    Q_soft(t+1,:) = Q_soft(t,:);
    
    % (2) update value estimate of chosen action: Q(new) = Q(old) + alpha * (R - Q(old))
    Q_soft(t+1,A_soft(t)) = Q_soft(t,A_soft(t)) + alpha * (R_soft(t) - Q_soft(t,A_soft(t)));
    
end

%% plot stuff

% % compute cumulative reward curves
cumulR_soft = cumsum(R_soft);
% % cumulR_epsilon = cumsum(R_epsilon);
% 
% % compute action frequency histograms
% [Nsoft,Esoft] = histcounts(A_soft,N,'BinMethod','integers');
% 
% subplot(2,1,1);plot(Q_soft)
% ylabel('Q: softmax')
% legend(labels)
% subplot(2,1,2)
% stairs(A_soft); set(gca,'YLim',[0.5 N+0.5])
% ylabel('Chosen action: softmax'); xlabel('Trials')
% 
% figure; stairs(cumulR_soft);
% xlabel('Trial'); ylabel('Accumulated reward')
% 
% figure;
% histogram(A_soft,N,'BinMethod','integers','Normalization','probability');
% set(gca,'YLim',[0 1])
% xlabel('Action'); ylabel('Proportion: softmax')

%% nested function
    function out = exp_beta(A,tau)
        x = 0:T;
        %tau = 10 % decay rate
        out = A*beta.^(-x/tau)
    end

    function out = lin_beta(rate)
        x = 0:T;
        out = -rate*x+beta
    end

end