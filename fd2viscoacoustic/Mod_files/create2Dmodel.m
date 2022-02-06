clear; clc; close all;

fp = fopen("../para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

nx    = A{1}(2);
nz    = A{1}(3);



c0   = 1800;
c1   = 2200;


c            = ones(nz,nx) * c0;
c(200:end,:) = c1;


fp = fopen("rho.bin","w");
fwrite(fp,c,"float32");
fclose(fp);



