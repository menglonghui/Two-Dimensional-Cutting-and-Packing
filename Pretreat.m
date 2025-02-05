function [Parent,template,I2_t,newRc] = Pretreat(I1,I2,nPop)
%PRETREAT 此处显示有关此函数的摘要
%   此处显示详细说明
% 生成空的结构体
template.x = [];
template.y = [];
template.z = [];
template.tu = [];
template.ddx = [];
template.ddy = [];
template.sita = [];

Parent=repmat(template, nPop, 1);% B = repmat(A,m,n)，将矩阵 A 复制 m×n 块
I2=trans2mid(I2,0);
newRc=I1;%I1最初为空白板
I2_t=I2;
%     sw=s_w(I2_t);%shortest width
end

