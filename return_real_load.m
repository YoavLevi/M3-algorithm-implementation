function [real_load] = return_real_load(q_i,real_q_i)

csum = cumsum(q_i);
csum = [0,csum];
real_load = zeros(1,numel(q_i));
for i=1:numel(q_i)
    real_load(i) = sum(real_q_i(csum(i)+1:csum(i+1)));
end
% return
end
