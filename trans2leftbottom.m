function Imid = trans2leftbottom(I)
%��ͼ��������м�λ��

    nl=0;%��߿ճ���
    nb=0;%�±߿ճ���

    for i=1:numel(I(1,:)) %newRec�ĳ�ʼֵΪI1
        if sum(I(:,i))==0
            nl=nl+1;
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

    Imid=imtranslate(I,[-nl, nb]);

end


