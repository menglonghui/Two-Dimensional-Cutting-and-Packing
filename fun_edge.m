function [I_final] = fun_edge(I)
%FUN_EDGE 此处显示有关此函数的摘要
%   此处显示详细说明
I_edge = edge(I,'Canny');

% curve = drawpolygon('LineWidth',1,'Color','r');

mask = imfill(I_edge,'holes');
figure,imshow(mask);


end

