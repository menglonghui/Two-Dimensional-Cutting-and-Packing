function [I_rec] = shape_to_rec(I)
%SHAPE 此处显示有关此函数的摘要
%   此处显示详细说明
lx=l_x(I);
ly=l_y(I);

[nl,~,nt,~] = margin(I);

for i = (nt + 1) : (nt + ly)
    for j = (nl + 1) : (nl + lx)        
        I(i,j) = 1;        
    end
end

I_rec = I;