function [Parent,template,I2_t,newRc] = Pretreat(I1,I2,nPop)
%PRETREAT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% ���ɿյĽṹ��
template.x = [];
template.y = [];
template.z = [];
template.tu = [];
template.ddx = [];
template.ddy = [];
template.sita = [];

Parent=repmat(template, nPop, 1);% B = repmat(A,m,n)�������� A ���� m��n ��
I2=trans2mid(I2,0);
newRc=I1;%I1���Ϊ�հװ�
I2_t=I2;
%     sw=s_w(I2_t);%shortest width
end

