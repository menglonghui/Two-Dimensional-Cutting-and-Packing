function Parent = Initial_pop_mix(Parent,nPop,nVar)
%INITIAL_POP �˴���ʾ�йش˺�����ժҪ
%������ɳ�ʼ��Ⱥ ���ü����ʼ��ȺҲ������������
for i = 1 : nPop    % nPop����Ⱥ��ģ��С
    
    Parent(i).x = randi([0, 1], 1, nVar); % nVar��x�ĳ���
    
end
end

