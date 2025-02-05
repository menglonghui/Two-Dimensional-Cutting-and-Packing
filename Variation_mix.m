function [Offspring] = Variation_mix(Offspring,heiban)
%VARIATION_MIX 此处显示有关此函数的摘要
%   此处显示详细说明
L_nR=numel(heiban(:,1))+numel(heiban(1,:));
for k = 1 : numel(Offspring) % 变异操作，nC为种群数量
    Offspring(k).x = mutatePop_mix(Offspring(k).x,L_nR);%x为二进制的位数  
end
end
