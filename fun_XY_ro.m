function Yx = fun_XY_ro(I,sita0)
%FUN_Y �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   �˴���ʾ��ϸ˵��
%     nl=0;%��߿ճ���
%     nr=0;%�ұ߿ճ���
%     nt=0;%�ϱ߿ճ���
%     nb=0;%�±߿ճ���
% NewR=zeros(round(0.8*l_w(I)));
% 
% for i=1:numel(I(:,1))
%     for j=1:numel(I(1,:))
%         NewR(i,j)=I(i,j);
%     end
% end
% 
% I=NewR;

I=trans2mid(I,0);
I_ro=imrotate(I,sita0,'bicubic');
% imshow(I_ro);
Yx=l_x(I_ro)*1000 + l_x(I)/10 + l_y(I)/10000;
%  Yx=l_x(I_ro);
end

