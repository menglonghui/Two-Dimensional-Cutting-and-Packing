function Yx = fun_YX(I)
%FUN_Y �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   �˴���ʾ��ϸ˵��
%     nl=0;%��߿ճ���
%     nr=0;%�ұ߿ճ���
%     nt=0;%�ϱ߿ճ���
%     nb=0;%�±߿ճ���
% 
%     for i=1:numel(I(1,:)) %newRec�ĳ�ʼֵΪI1
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

