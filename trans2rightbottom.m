function Imid = trans2rightbottom(I)
%��ͼ��������м�λ��

    nr=0;%��߿ճ���
    nb=0;%�±߿ճ���

    for i=numel(I(1,:)):-1:1 %newRec�ĳ�ʼֵΪI1
        if sum(I(:,i))==0
            nr=nr+1;
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

    Imid=imtranslate(I,[nr, nb]);

end


