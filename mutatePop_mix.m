function p = mutatePop_mix(x,L_nR)  %����
%MUTATEPOP_MIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
max_end=max(2,24-floor(L_nR/220));
min=1;
max_value=randi([1,max_end]);
byws=randi([min, max_value])-randi([0,1]);%����λ��
m = numel(x);
for j=1:byws
    s = randi([1, m]);
    x(s)=1-x(s);
end
p = x;
end
