function Parent = Initial_pop(Parent,nPop,nVar,tx)
%INITIAL_POP 此处显示有关此函数的摘要
%随机生成初始种群 被裁剪后初始种群也是在这里生成
for i = 1 : nPop    % nPop，种群规模大小
    
    Parent(i).x = randi([0, 1], 1, nVar); % nVar：x的长度
    if tx == 0
        Parent(i).x(end) = 1 - (Parent(i).x(25));
    else
        Parent(i).x(end) = (Parent(i).x(25));
    end
    
end
end

