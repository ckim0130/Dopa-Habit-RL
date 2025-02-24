function [time_1, habit_num] = bandit_two(n_ban,T,alpha,beta)
%% Agent parameters
% [temp,max_arm,j,k] = bandit(n_ban,T,alpha,beta)
% bandit(4,1000,0.1,2)
close all;

%% bandit parameters
% key is expected reward E(i) = P(i) x R(i)
% P = [1 1 1 1];              % deterministic bandit
% R = [0.4 0.7 0.1 0.5];      % variable rewards

P = ones(1,n_ban);
% R = int32(randi([0,1],[1,n_ban])) Question: 0 or 1??
% R = rand(1,n_ban)
R = [0.2 0.7];
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


%% dynamic beta
% dynamic_beta = exp_beta(1.3,180);%exp_beta(A,tau) % t = 1000
% figure; plot(dynamic_beta)
% dynamic_beta = exp_beta(1.3,18); % t = 100
% dynamic_beta = sig_beta(100,0.1);
dynamic_beta = heavi_beta(80,500);
% dynamic_beta = relu_beta(0.2);
% figure; plot(dynamic_beta)
% dynamic_beta = lin_beta(1.4);

% deviation = 0.2; %goal-directed 
deviation = 0.01; % small deviation
temp=[];
count = 0;
time = [];

for t = 1:T
    % choose action (pull arm of bandit) based on current estimated value of each bandit
    %     P_action_soft = exp(Q_soft(t,:)*dynamic_exp(t)) ./ sum(exp(Q_soft(t,:)*dynamic_exp(t))); % softmax P(action)
    P_action_soft = exp(Q_soft(t,:)*dynamic_beta(t)) ./ sum(exp(Q_soft(t,:) * dynamic_beta(t))); % softmax P(action)
    A_soft(t) = find(rand <= cumsum(P_action_soft),1);  % choose action based on those probabilities
    
    % was action rewarded?
    R_soft(t) = (rand <= P(A_soft(t))) * R(A_soft(t)); 
    % Update values of actions giiven reward:
    % (1) all actions not chosen keep same value estimates
    Q_soft(t+1,:) = Q_soft(t,:);
    
    % (2) update value estimate of chosen action: Q(new) = Q(old) + alpha * (R - Q(old))
    Q_soft(t+1,A_soft(t)) = Q_soft(t,A_soft(t)) + alpha * (R_soft(t) - Q_soft(t,A_soft(t)));

    if Q_soft(t,2) > R(2) - deviation
        R(2) = 0.1;
        time = [time, t];
%         count = count + 1;
%         if count == 700
%             R(2) = 0.1;
%             time = [time,t];
%         end
    end
    maxx = Q_soft(t,:); % max arm for a given step
    max_arm = find(ismember(maxx, max(maxx)));
    temp = [temp,max_arm];

    
end

time_1 = time(1);
% testing how many it will continut to press until the next one
habit = temp(min(time):end);

for habit_num = 1:length(habit)
    if habit(habit_num) ~= 2
        break
    end
end

habit_num

% time interval of 50
% habit = temp(min(time):min(time)+50);

% for c = 1:N
%     one = test(test(:,1)==c,:)
%     two = test(test(:,1)==c,:)
%     three = test(test(:,1)==c,:)
%     four = test(test(:,1)==c,:)
% end
% 
% asymp_1 = one(1,2);
% asymp_2 = two(1,2);
% asymp_3 = three(1,2);
% asymp_4 = four(1,2);

%     maxx = Q_soft(t,:); % max arm for a given step
%     max_arm = find(ismember(maxx, max(maxx)));
%     temp = [temp,max_arm];


%% when the highest arm surpasses the other arms

% arm_temp = [];
% for c = 1:N
%     arm_temp = [arm_temp,sum(A_soft == c)];
%     max_arm = find(ismember(arm_temp, max(arm_temp(:))));
% end
% 
% 
% % highest reward among n_bandit arms
% for j = 2:T % starting from 2 bc 1 is noisy
%     Q_soft_row = Q_soft(j,:);
%     if Q_soft_row(max_arm) == max(Q_soft_row) 
%         sprintf('highest reward: %d', j);
%         break
%     end
% end
% 
% % asymptote
% deviation = 0.01;
% 
% for k = 2:T
%     if Q_soft(k,max_arm) > R(max_arm) - deviation % reach asymptote
%         sprintf('Asymptote: %d', k);
%         break
%     end
% end

% high_arm_val = find(ismember(Q_soft(:,max_arm), max(Q_soft(:,max_arm))));
% min_high_arm = min(high_arm_val)


    
% %% plot stuff
% 
% % compute cumulative reward curves
% cumulR_soft = cumsum(R_soft);
% cumulR_max = cumulR_soft(end);


% 
% % compute action frequency histograms
% [Nsoft,Esoft] = histcounts(A_soft,N,'BinMethod','integers');
% 
% subplot(2,1,1);plot(Q_soft);
% ylabel('Q: softmax');
% legend(labels);
% subplot(2,1,2);
% stairs(A_soft); set(gca,'YLim',[0.5 N+0.5]);
% ylabel('Chosen action: softmax'); xlabel('Trials');
% 
% figure; stairs(cumulR_soft);
% xlabel('Trial'); ylabel('Accumulated reward');
% 
% 
% figure;
% % histogram(A_soft,N,'Normalization','probability');
% histogram(A_soft,N,'BinMethod','integers','Normalization','probability')
% % histogram(A_soft);
% set(gca,'YLim',[0 1]);
% xlabel('Action'); ylabel('Proportion: softmax');

%% nested function
    function out = exp_beta(A,tau)
        x = 0:T;
        %tau = 10 % decay rate
        out = A*beta.^(x/tau);
    end

    function out = lin_beta(rate)
        x = 0:T;
        out = rate*x+beta;
    end

    function out = sig_beta(A,b)
        x = -T/2:T/2;
        out = A./(1+exp(-b*x))+2;
    end

    function out = heavi_beta(B,onset) %heaviside function  
        x = 0:T;
        for i = 1:length(x)
            if x(i) > onset
                out(i) = B;
            else
                out(i) = 2;
            end
        end      
    end

    function out = relu_beta(rate)
        x = -T/2:T/2;
        out_temp = rate*x+beta;
        out = max(2,out_temp);
    end
end