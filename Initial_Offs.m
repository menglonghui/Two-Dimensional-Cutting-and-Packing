function Offspring = Initial_Offs(template,Parent,nC)
%INITIAL_OFFS �˴���ʾ�йش˺�����ժҪ
    Offspring = repmat(template, nC/2, 2); %�����Ӵ��ṹ��
    for j = 1 : nC / 2     % �������
        p1 = selectPop(Parent);% ѡ����������ѡȡ������ѡ������
        p2 = selectPop(Parent);% ѡ����������ѡȡ������ѡ������
        [Offspring(j, 1).x, Offspring(j, 2).x] = crossPop(p1.x, p2.x);%��ѡ���ĸ�����н��棬Ȼ��ֱ���Ϊ�Ӵ�
    end
    Offspring = Offspring(:);  %��Offspringת����һά����
end

