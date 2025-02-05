function S = fun_S(I)
%FUN_S 此处显示有关此函数的摘要
%   此处显示详细说明
%     nl=0;%左边空出的
%     nr=0;%右边空出的
%     nt=0;%上边空出的
%     nb=0;%下边空出的

%     I_t=imrotate(I,sita,'crop');
%     shape=zeros(180,1);

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
%     S=(numel(I(1,:))-nl-nr)*(numel(I(:,1))-nt-nb);
S=l_x(I)*l_y(I);
end

