clear; clc; close all;

nx     = 500;
nz     = 500;
gWidth = 100;

nx     = nx +  2 * gWidth;

c0   = 2000;
c1   = 1500;
cutU = 300;
cutD = 400;

c  = ones(nz,nx) * c0;
c(cutU:cutD,:) = c1;

fp = fopen("c.bin","w");
fwrite(fp,c,"float32");
fclose(fp);



