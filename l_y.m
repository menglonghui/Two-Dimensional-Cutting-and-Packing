function ly = l_y(I)
%L_W �˴���ʾ�йش˺�����ժҪ
    nt=0;
    nb=0;

    for i=1:numel(I(:,1)) %newRec�ĳ�ʼֵΪI1
        if sum(I(i,:))==0
            nt=nt+1;
        else
            break
        end
    end  %���ڷ���ɹ����ͼ�񣬼���nl��nr��nt��nb�ĸ�ֵ

    for i=numel(I(:,1)):-1:1
        if sum(I(i,:))==0
            nb=nb+1;
        else
            break
        end
    end
    
    ly=numel(I(:,1))-nb-nt;
    
end

