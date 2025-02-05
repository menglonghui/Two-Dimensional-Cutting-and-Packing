function [newRc,I2_t,I_last,Parent_last,nl_last,nr_last,nb_last,nt_last] = Cropping(I,newRc,I2_t,Parent)
%CROPPING 此处显示有关此函数的摘要
%   此处显示详细说明
    nl=0;%左边空出的
    nr=0;%右边空出的
    nt=0;%上边空出的
    nb=0;%下边空出的

    %ita=Nw_imerged/(sum(sum(I))+sum(sum(~I)));

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

    for i=1:numel(I(:,1))
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

    %下面开始重新定义矩形
    newRc=zeros(numel(newRc(:,1))-nt-nb,numel(newRc(1,:))-nl-nr);
    %newRc为新创建的一个与成功分离后的图像经过裁剪后大小一致的空白板

    %[I2_t_ny, I2_t_nx]=size(newRc);
    [I_ny, I_nx]=size(I);
    I2_t=imtranslate(I2_t,[round((nl-nr)/2),round((nt-nb)/2)]);   %I2_t_nx表示横向像素的数目
    %裁剪前先对I2_t图案的位置进行调整，确保裁剪后图案还是在中间的位置
    I2_t_temp=zeros(numel(newRc(:,1)),numel(newRc(1,:)));

    for i=1:I_ny-nt-nb
        for j=1:I_nx-nl-nr
            I2_t_temp(i,j)=I2_t(nt+i,nl+j);
        end
    end

    I2_t=I2_t_temp;

    I2_t=trans2mid(I2_t);

    %因为分离成功了，所以下面保存该次的搜索结果，以防下次无法找到合适的解，就返回这次的结果
    I_last=I;
    Parent_last=Parent;
    nl_last=nl;
    nr_last=nr;
    nb_last=nb;
    nt_last=nt;


end


