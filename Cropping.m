function [newRc,I2_t,I_last,Parent_last,nl_last,nr_last,nb_last,nt_last] = Cropping(I,newRc,I2_t,Parent)
%CROPPING �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    nl=0;%��߿ճ���
    nr=0;%�ұ߿ճ���
    nt=0;%�ϱ߿ճ���
    nb=0;%�±߿ճ���

    %ita=Nw_imerged/(sum(sum(I))+sum(sum(~I)));

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

    %���濪ʼ���¶������
    newRc=zeros(numel(newRc(:,1))-nt-nb,numel(newRc(1,:))-nl-nr);
    %newRcΪ�´�����һ����ɹ�������ͼ�񾭹��ü����Сһ�µĿհװ�

    %[I2_t_ny, I2_t_nx]=size(newRc);
    [I_ny, I_nx]=size(I);
    I2_t=imtranslate(I2_t,[round((nl-nr)/2),round((nt-nb)/2)]);   %I2_t_nx��ʾ�������ص���Ŀ
    %�ü�ǰ�ȶ�I2_tͼ����λ�ý��е�����ȷ���ü���ͼ���������м��λ��
    I2_t_temp=zeros(numel(newRc(:,1)),numel(newRc(1,:)));

    for i=1:I_ny-nt-nb
        for j=1:I_nx-nl-nr
            I2_t_temp(i,j)=I2_t(nt+i,nl+j);
        end
    end

    I2_t=I2_t_temp;

    I2_t=trans2mid(I2_t);

    %��Ϊ����ɹ��ˣ��������汣��ôε�����������Է��´��޷��ҵ����ʵĽ⣬�ͷ�����εĽ��
    I_last=I;
    Parent_last=Parent;
    nl_last=nl;
    nr_last=nr;
    nb_last=nb;
    nt_last=nt;


end


