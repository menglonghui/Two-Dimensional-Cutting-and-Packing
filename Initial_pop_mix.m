function Parent = Initial_pop_mix(Parent,nPop,nVar)
%INITIAL_POP 此处显示有关此函数的摘要
%随机生成初始种群 被裁剪后初始种群也是在这里生成
for i = 1 : nPop    % nPop，种群规模大小
    
    Parent(i).x = randi([0, 1], 1, nVar); % nVar：x的长度
    
end
end

