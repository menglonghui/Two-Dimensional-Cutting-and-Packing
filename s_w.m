function sw = s_w(I)
%WIDTH_SHORTEST 此处显示有关此函数的摘要
sw=numel(I(1,:));
    for sita=0:90:180
        nl=0;%左边空出的
        nr=0;%右边空出的
        nt=0;%上边空出的
        nb=0;%下边空出的
        I=imrotate(I,sita,'crop','bicubic');
        for i=1:numel(I(1,:)) %newRec的初始值为I1
            if sum(I(:,i))==0
                nl=nl+1;
            else
                break
            end
        end  %基于分离成功后的图像，计算nl，nr，nt，nb四个值

        for i=numel(I(1,:)):-1:1
            if sum(I(:,i))==0
                nr=nr+1;
            else
                break
            end
        end

        for i=1:numel(I(:,1))%numel(I(:,1))为图像的上下距离（像素）
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

