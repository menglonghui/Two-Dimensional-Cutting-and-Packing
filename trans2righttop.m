function Imid = trans2righttop(I)
%��ͼ��������м�λ��

    nr=0;%��߿ճ���
    nt=0;%�±߿ճ���

    for i=numel(I(1,:)):-1:1 %newRec�ĳ�ʼֵΪI1
        if sum(I(:,i))==0
            nr=nr+1;
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

    Imid=imtranslate(I,[nr, -nt]);

end


