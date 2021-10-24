clear; clc; close all;


fp = fopen("para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

%Parsing parameters
nt    = A{1}(1);
nx    = A{1}(2);
nz    = A{1}(3);
isnap = A{1}(9);
gWidth= A{1}(13);

nx    = nx + 2 * gWidth;
ns    = round(nt / isnap);
f = fopen('OUTPUT/field_1.bin','r');
U = fread(f, 'float64');
fclose(f);

U = reshape(U,ns,nz,nx);


maxamp = max(max(max(U)));  % max amplitude calculation for clip
clip = 96;  % clip value
ampclip = (1-clip/100)*maxamp;
data    = U; 
data(data > ampclip)= ampclip; % clipping positive
data(data < -ampclip)= -ampclip; % clipping positive


for is=1:ns
    imagesc(reshape(data(is,:,:),nz,nx));   
    colormap(redblue);
    pause(.099);
end

