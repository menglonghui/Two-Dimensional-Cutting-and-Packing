function new = selectPop(Parent)  %ѡ�������Ӹ�����Ⱥ�н���ѡ��
%���̶�ѡ��
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


%������ѡ��
%     n = numel(Parent);
%     index = randperm(n);  %y = randperm(n)    y�ǰ�1��n��Щ��������ҵõ���һ���������С�
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