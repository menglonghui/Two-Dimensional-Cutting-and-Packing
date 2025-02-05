function Yx = fun_YX(I)
%FUN_Y 此处显示有关此函数的摘要
%   此处显示详细说明
%   此处显示详细说明
%     nl=0;%左边空出的
%     nr=0;%右边空出的
%     nt=0;%上边空出的
%     nb=0;%下边空出的
% 
%     for i=1:numel(I(1,:)) %newRec的初始值为I1
%         if sum(I(:,i))==0
%             nl=nl+1;
%         else
%             break
%         end
%     end
% 
%     for i=numel(I(1,:)):-1:1
%         if sum(I(:,i))==0
%             nr=nr+1;
%         else
%             break
%         end
%     end
% 
%     for i=1:numel(I(:,1))
%         if sum(I(i,:))==0
%             nt=nt+1;
%         else
%             break
%         end
%     end
% 
%     for i=numel(I(:,1)):-1:1
%         if sum(I(i,:))==0
%             nb=nb+1;
%         else
%             break
%         end
%     end
    
%     Yx=(numel(I(:,1))-nt-nb)*5+(numel(I(1,:))-nl-nr)/5;
%     Yx=(numel(I(:,1))-nt-nb)/5+(numel(I(1,:))-nl-nr)*5;
% Yx=(numel(I(:,1))-nt-nb)+(numel(I(1,:))-nl-nr);
    Yx=l_y(I)*5+l_x(I)/5;
end

