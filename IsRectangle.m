function [isrectangle] = IsRectangle(I)
%ISRECTANGLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I = imfill(I,'holes');
lx = l_x(I);
ly = l_y(I);
N_p = sum(I(:));

if (lx * ly > N_p * 0.95) && (lx * ly < N_p * 1.05)
    isrectangle = 1;
else
    isrectangle = 0;
end

end

