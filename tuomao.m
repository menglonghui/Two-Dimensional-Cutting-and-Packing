function I_final = tuomao(I)
%TUOMAO 此处显示有关此函数的摘要
%   此处显示详细说明
% I_final=zeros(numel(I(:,1)),numel(I(1,:)));
% I_final=I;

I_base_o = zeros(size(I));
I_base_expand = zeros(2*size(I));
I = trans2lefttop(I);
I = imcrop(I,[0,0,l_x(I),l_y(I)]);
I = paste_img(I_base_expand,I,1,1);
I = trans2mid(I,0);
I_temp=I;
lx=l_x(I);
ly=l_y(I);
leftpoint=round((numel(I(1,:))-lx-2)/2);
rightpoint=numel(I(1,:))-leftpoint;
toppoint=round((numel(I(:,1))-ly-2)/2);
bottompoint=numel(I(:,1))-toppoint;
for i=toppoint:bottompoint
    for j=leftpoint:rightpoint
        %下面进行脱毛
        if I(i,j)==1 && (I(i-1,j-1) + I(i-1,j) + I(i-1,j+1) + I(i,j-1) + I(i,j+1) + I(i+1,j-1) + I(i+1,j) + I(i+1,j+1) <= 3) &&...
                ((I(i-1,j)==0 && I(i+1,j)==0) || (I(i,j-1)==0 && I(i,j+1)==0))
            I_temp(i,j)=0;
        end
    end
end

I = I_temp;
I = trans2lefttop(I);
I = imcrop(I,[0,0,l_x(I),l_y(I)]);
I = paste_img(I_base_o,I,1,1);
I = trans2mid(I,0);
I_final = I;

end
