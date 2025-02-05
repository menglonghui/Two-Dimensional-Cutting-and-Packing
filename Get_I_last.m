function I = Get_I_last(I_last,nl_last,nr_last,nb_last,nt_last)
%GET_I_LAST 此处显示有关此函数的摘要
%   此处显示详细说明
    [ny,nx]=size(I_last);
    I=zeros(ny-nb_last-nt_last,nx-nl_last-nr_last);

    for i=1:(ny-nb_last-nt_last)
        for j=1:(nx-nl_last-nr_last)
            I(i,j)=I_last(i+nt_last,j+nl_last);
        end
    end

end

