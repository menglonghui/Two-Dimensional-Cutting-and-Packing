function [dx, dy] = pop_decode_xy(bi,Nx,Ny,lx,ly)
%�ú����������Ʊ���ת��Ϊʮ���Ʊ���
    dx=(Nx-lx)*(bin2dec(num2str(bi(:,1:12)))/4095);%��ΧΪ0~(Nx-sw)
    dy=(Ny-ly)*(bin2dec(num2str(bi(:,13:24)))/4095);%��ΧΪ0~(Ny-sw)
    %     sita=bin2dec(num2str(bi(:,25)));

    dx=round((Nx-lx)/2-dx);%dx��ȡֵ��Χ
    dy=round((Ny-ly)/2-dy);%dy��ȡֵ��Χ
    %     sita=sita*180;
end