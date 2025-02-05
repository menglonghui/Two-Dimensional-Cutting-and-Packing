function [Offspring] = Variation(Offspring,~,nn,newRc,tx)
%VARIATION  此处显示有关此函数的摘要
%变异操作（同时计算所有的元素的y值）
%tx 是否同向
    L_nR=numel(newRc(:,1))+numel(newRc(1,:));
    for k = 1 : numel(Offspring) % 变异操作，nC为种群数量
        Offspring(k).x = mutatePop(Offspring(k).x,nn,L_nR);%x为二进制的位数
        if tx == 0
            Offspring(k).x(50) = 1 - (Offspring(k).x(25));%保持两者角度相差180度
        else
            Offspring(k).x(50) = Offspring(k).x(25);
        end
    end
end