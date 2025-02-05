function new = selectPop(Parent)  %选择函数，从父代种群中进行选择
%轮盘赌选择
flag=0;
popsize = numel(Parent);%
% i=1;
fitness=zeros(popsize,1);
for k=1:popsize
    fitness(k)=Parent(k).y;
end
p=fitness/sum(fitness);
Cs=cumsum(p);
%R=sort(rand(popsize,1));


i=1;
while i<=popsize
    if rand<=Cs(i)
        new=Parent(i);
        flag=1;
        break
    else
        i=i+1;
    end
end

if flag==0
    i=fix(rand*popsize)+1;
    new=Parent(i);
end


%锦标赛选择法
%     n = numel(Parent);
%     index = randperm(n);  %y = randperm(n)    y是把1到n这些数随机打乱得到的一个数字序列。
%     p1 = Parent(index(1));
%     p2 = Parent(index(2));
%
%     if p1.y >= p2.y
%         p = p1;
%     else
%         p = p2;
%     end
%
%     new=Parent(p);

end