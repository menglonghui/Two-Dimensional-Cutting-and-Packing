function Imid = trans2lefttop(I)
%��ͼ��������м�λ��

    nl=0;%��߿ճ���
    nt=0;%�±߿ճ���

    for i=1:numel(I(1,:)) %newRec�ĳ�ʼֵΪI1
        if sum(I(:,i))==0
            nl=nl+1;
        else
            break
        end
    end

    for i=1:numel(I(:,1))
        if sum(I(i,:))==0
            nt=nt+1;
        else
            break
        end
    end

    Imid=imtranslate(I,[-nl, -nt]);

end


