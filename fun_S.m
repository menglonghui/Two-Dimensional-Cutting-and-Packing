function S = fun_S(I)
%FUN_S �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%     nl=0;%��߿ճ���
%     nr=0;%�ұ߿ճ���
%     nt=0;%�ϱ߿ճ���
%     nb=0;%�±߿ճ���

%     I_t=imrotate(I,sita,'crop');
%     shape=zeros(180,1);

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
%     S=(numel(I(1,:))-nl-nr)*(numel(I(:,1))-nt-nb);
S=l_x(I)*l_y(I);
end

