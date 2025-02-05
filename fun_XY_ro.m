function Yx = fun_XY_ro(I,sita0)
%FUN_Y 此处显示有关此函数的摘要
%   此处显示详细说明
%   此处显示详细说明
%     nl=0;%左边空出的
%     nr=0;%右边空出的
%     nt=0;%上边空出的
%     nb=0;%下边空出的
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

