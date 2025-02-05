function Imid = trans2leftbottom(I)
%将图像调整至中间位置

    nl=0;%左边空出的
    nb=0;%下边空出的

    for i=1:numel(I(1,:)) %newRec的初始值为I1
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


