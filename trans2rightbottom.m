function Imid = trans2rightbottom(I)
%将图像调整至中间位置

    nr=0;%左边空出的
    nb=0;%下边空出的

    for i=numel(I(1,:)):-1:1 %newRec的初始值为I1
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


