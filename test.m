tt = lin_beta(2);
t1 = max(2,tt);
figure; plot(tt);figure; plot(t1)

function out = lin_beta(rate)
        x = -50:50;
        out = rate*x+2;
end
    

% function out = relu_beta(x,beta)
%         x = 0:100;
%         for i = 1:length(x)
%             if x(i) < 50
%                 out(i) = 0;
%             else
%                 out(i) = 1*x + beta;
%             end
%         end
% end
   