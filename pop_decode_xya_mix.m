function [dx, dy, sita] = pop_decode_xya_mix(bi,Nx,Ny,lx,ly)
%�ú����������Ʊ���ת��Ϊʮ���Ʊ���
    dx = (Nx - lx) * (bin2dec(num2str(bi(01 : 12))) / 4095);%��ΧΪ0~(Nx-sw)
    dy = (Ny - ly) * (bin2dec(num2str(bi(15 : 26))) / 4095);%��ΧΪ0~(Ny-sw)
    sita = bin2dec(num2str(bi(13:14)));

    dx = round(1.01 * ((Nx - lx) / 2 - dx));%dx��ȡֵ��Χ
    dy = round(1.01 * ((Ny - ly) / 2 - dy));%dy��ȡֵ��Χ
    sita = sita * 90;
end
