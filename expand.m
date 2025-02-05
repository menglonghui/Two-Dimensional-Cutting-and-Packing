function [I_final] = expand(I,partgap)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
% I_base = zeros(size(I));
% I = trans2lefttop(I);
% I = imcrop(I,[0,0,l_x(I),l_y(I)]);
%
% % I_base = zeros(l_y(I)+2*partgap,l_x(I)+2*partgap);
% I = paste_img(I_base,I,1,1);
% I = trans2mid(I,0);
I_base = zeros(numel(I(:,1)) + partgap, numel(I(1,:)) + partgap);
% I_base = zeros(numel(I(:,1)), numel(I(1,:)));
I = paste_img(I_base,I,1,1);
I = trans2mid(I,0);
I_final = I_base;
% I_final = zeros(size(I));
lx = l_x(I);
ly = l_y(I);


leftpoint = round((numel(I(1,:))-lx-4)/2);
rightpoint = numel(I(1,:))-leftpoint;
toppoint = round((numel(I(:,1))-ly-4)/2);
bottompoint = numel(I(:,1))-toppoint;
% imshow(I)
for i = toppoint : bottompoint
    for j = leftpoint : rightpoint
        N_around = I(i-1,j-1) + I(i-1,j) + I(i-1,j+1) + I(i,j-1) + I(i,j+1) + I(i+1,j-1) + I(i+1,j) + I(i+1,j+1);
        
        if I(i,j) == 1 && (N_around < 3)
            
            for p = (i - ceil(partgap / 2)) : (i + ceil(partgap / 2))
                for q = (j - ceil(partgap / 2)) : (j + ceil(partgap / 2))
                    %                     if ((p - i)^2 + (q - j)^2) <= ((partgap / 2)^2)
                    if ((p - i)^2 + (q - j)^2) <= ((partgap / 2 + 1)^2)
                        I_final(p,q) = 1;
                    end
                end
            end
            
        elseif I(i,j) == 1 && (N_around >= 3 && N_around < 8)
            
            for p = (i - ceil(partgap / 2)) : (i + ceil(partgap / 2))
                for q = (j - ceil(partgap / 2)) : (j + ceil(partgap / 2))
                    %                     if ((p - i)^2 + (q - j)^2) <= ((partgap / 2)^2)
                    if ((p - i)^2 + (q - j)^2) <= ((partgap / 2 + 1)^2 - 1)
                        I_final(p,q) = 1;
                    end
                end
            end
        end
        
    end
end
I_final = I_final | I;
% imshow(I_final)
end
