function [Offspring] = Variation_mix(Offspring,heiban)
%VARIATION_MIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
L_nR=numel(heiban(:,1))+numel(heiban(1,:));
for k = 1 : numel(Offspring) % ���������nCΪ��Ⱥ����
    Offspring(k).x = mutatePop_mix(Offspring(k).x,L_nR);%xΪ�����Ƶ�λ��  
end
end
