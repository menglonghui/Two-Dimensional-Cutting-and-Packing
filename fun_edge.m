function [I_final] = fun_edge(I)
%FUN_EDGE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I_edge = edge(I,'Canny');

% curve = drawpolygon('LineWidth',1,'Color','r');

mask = imfill(I_edge,'holes');
figure,imshow(mask);


end

