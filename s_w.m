function sw = s_w(I)
%WIDTH_SHORTEST �˴���ʾ�йش˺�����ժҪ
sw=numel(I(1,:));
    for sita=0:90:180
        nl=0;%��߿ճ���
        nr=0;%�ұ߿ճ���
        nt=0;%�ϱ߿ճ���
        nb=0;%�±߿ճ���
        I=imrotate(I,sita,'crop','bicubic');
        for i=1:numel(I(1,:)) %newRec�ĳ�ʼֵΪI1
            if sum(I(:,i))==0
                nl=nl+1;
            else
                break
            end
        end  %���ڷ���ɹ����ͼ�񣬼���nl��nr��nt��nb�ĸ�ֵ

        for i=numel(I(1,:)):-1:1
            if sum(I(:,i))==0
                nr=nr+1;
            else
                break
            end
        end

        for i=1:numel(I(:,1))%numel(I(:,1))Ϊͼ������¾��루���أ�
            if sum(I(i,:))==0
                nt=nt+1;
            else
                break
            end
        end

        for i=numel(I(:,1)):-1:1
            if sum(I(i,:))==0
                nb=nb+1;
            else
                break
            end
        end
        
        ws_temp=min(numel(I(:,1))-nt-nb,numel(I(1,:))-nl-nr);
        sw=min(ws_temp,sw);

    end
end

