function [Offspring] = Variation(Offspring,~,nn,newRc,tx)
%VARIATION  �˴���ʾ�йش˺�����ժҪ
%���������ͬʱ�������е�Ԫ�ص�yֵ��
%tx �Ƿ�ͬ��
    L_nR=numel(newRc(:,1))+numel(newRc(1,:));
    for k = 1 : numel(Offspring) % ���������nCΪ��Ⱥ����
        Offspring(k).x = mutatePop(Offspring(k).x,nn,L_nR);%xΪ�����Ƶ�λ��
        if tx == 0
            Offspring(k).x(50) = 1 - (Offspring(k).x(25));%�������߽Ƕ����180��
        else
            Offspring(k).x(50) = Offspring(k).x(25);
        end
    end
end