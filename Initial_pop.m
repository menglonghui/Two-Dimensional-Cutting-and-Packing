function Parent = Initial_pop(Parent,nPop,nVar,tx)
%INITIAL_POP �˴���ʾ�йش˺�����ժҪ
%������ɳ�ʼ��Ⱥ ���ü����ʼ��ȺҲ������������
for i = 1 : nPop    % nPop����Ⱥ��ģ��С
    
    Parent(i).x = randi([0, 1], 1, nVar); % nVar��x�ĳ���
    if tx == 0
        Parent(i).x(end) = 1 - (Parent(i).x(25));
    else
        Parent(i).x(end) = (Parent(i).x(25));
    end
    
end
end

