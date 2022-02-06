clear; clc; close all;

fp = fopen("para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

nx    = A{1}(2);
nz    = A{1}(3);



c0   = 1;
c1   = 60;
c2   = 1000;


c  = ones(nz,nx) * c0;
%c(71:140,:) = c1;
%c(141:end,:) = c2;


fp = fopen("rho.bin","w");
fwrite(fp,c,"float32");
fclose(fp);



