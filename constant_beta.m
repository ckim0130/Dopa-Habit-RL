%%constant beta

temp = [];
% 
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
% %     maxx = Q_soft(t,:);
% %     max_arm = find(ismember(maxx, max(maxx)));
% %     temp = [temp,max_arm];
%     
% %     for c = 1:N
% %         test = [test,sum(A_soft == c)];
% %         max_arm = find(ismember(test, max(test(:))));
% %     end
% % 
% % 
% %     % highest reward among n_bandit arms
% %     for j = 2:T % starting from 2 bc 1 is noisy
% %         Q_soft_row = Q_soft(j,:);
% %         if Q_soft_row(max_arm) == max(Q_soft_row) 
% %             sprintf('highest reward: %d', j);
% %             break
% %         end
% %     end
%     
% end