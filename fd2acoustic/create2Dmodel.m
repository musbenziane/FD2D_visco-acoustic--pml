clear; clc; close all;

fp = fopen("para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

nx    = A{1}(2);
nz    = A{1}(3);



c0   = 2000;
c1   = 2000;
cutU = 600;
cutD = 1000;

c  = ones(nz,nx) * c0;
c(cutU:cutD,:) = c1;

fp = fopen("ctrial.bin","w");
fwrite(fp,c,"float32");
fclose(fp);



