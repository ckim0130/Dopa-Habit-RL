
T = 100
% results = zeros(T,2); % highest arm, highest value
check = zeros(T,T); %temp
% cumR = zeros(T,1);
for i = 1:T
%     [results(i,1), results(i,2)] = bandit_two(2,1000,0.1,2);
    check(i,:) = bandit(4,100,0.1,2);
%       cumR(i) = bandit(4,1000,0.1,20);
end
figure; surf(check); xlabel('Trials','FontSize',16); ylabel('Iterations','FontSize',16);set(gca,'FontSize',16);yticks([0 20 40 60 80 100]);yticklabels({'0','20','40','60','80','100'})
mean(results(:,2))
std(results(:,2))

x1 = sum(check(:) == 1);
x2 = sum(check(:) == 2);
x3 = sum(check(:) == 3);
x4 = sum(check(:) == 4);
final = x1+x3+x4

% reshape matrix into row
% row_vector = reshape(check.',1,[]); 

% arm1 = [];
% arm2 = [];
% arm3 = [];
% arm4 = [];
% 
% for j=1:T
%     k = results(j);
%     if k == 1
%         arm1 = [arm1;results(j,:)];
%     end
%     if k == 2
%         arm2 = [arm2;results(j,:)];
%     end
%     if k == 3
%         arm3 = [arm3;results(j,:)];
%     end
%     if k == 4
%         arm4 = [arm4;results(j,:)];
%     end
% end
% 
% %plot asymptote
% fig = figure;
% subplot(2,2,1); scatter(1:length(arm1),arm1(:,3)); title('arm1'); ylim([0 600])
% subplot(2,2,2); scatter(1:length(arm2),arm2(:,3)); title('arm2'); ylim([0 600])
% subplot(2,2,3); scatter(1:length(arm3),arm3(:,3)); title('arm3'); ylim([0 600])
% subplot(2,2,4); scatter(1:length(arm4),arm4(:,3)); title('arm4'); ylim([0 600])
% han = axes(fig,'visible','off');
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'asymptomatic value');
% xlabel(han,'# of times chosen');
% 
% 
% fig = figure;
% subplot(2,2,1); histogram(arm1(:,3)); title('arm1')
% subplot(2,2,2); histogram(arm2(:,3)); title('arm2') 
% subplot(2,2,3); histogram(arm3(:,3)); title('arm3')
% subplot(2,2,4); histogram(arm4(:,3)); title('arm4')
% han = axes(fig,'visible','off');
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'frequency');
% xlabel(han,'asymptomatic values');


% figure; stairs(cumulR_1); hold on;
% stairs(cumulR_2); hold on;
% legend('exponential','linear')