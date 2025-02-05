function Offspring = Initial_Offs(template,Parent,nC)
%INITIAL_OFFS 此处显示有关此函数的摘要
    Offspring = repmat(template, nC/2, 2); %创建子代结构体
    for j = 1 : nC / 2     % 交叉操作
        p1 = selectPop(Parent);% 选择操作，随机选取两个，选择优者
        p2 = selectPop(Parent);% 选择操作，随机选取两个，选择优者
        [Offspring(j, 1).x, Offspring(j, 2).x] = crossPop(p1.x, p2.x);%将选择后的个体进行交叉，然后分别作为子代
    end
    Offspring = Offspring(:);  %将Offspring转换成一维向量
end

