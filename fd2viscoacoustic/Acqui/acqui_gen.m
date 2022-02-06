clear; clc; close all;

nrcv  = 199;
drcv  = 6.25;
z0    = 12.5;

fp = fopen("../para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

%Parsing parameters
nx    = A{1}(2);
nz    = A{1}(3);
dx    = A{1}(4);
npml  = A{1}(7);


xrcv    = (npml-1)*dx:drcv:(nx+npml-1)*dx;
zrcv    = ones(length(xrcv),1) * z0;

f = fopen("rcv.acqui", 'w');

for i = 1:length(xrcv)
fprintf(f, '%.2f %.2f \n',zrcv(i), xrcv(i));
end

fclose(f);