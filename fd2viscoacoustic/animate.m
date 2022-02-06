clear; clc; close all;


fp = fopen("para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

%Parsing parameters
nt    = A{1}(1);
nx    = A{1}(2);
nz    = A{1}(3);
npml  = A{1}(7);
isnap = A{1}(10);
ns    = round(nt / isnap);

nx = nx + npml*2;
nz = nz + npml;

f = fopen('OUTPUT/field_P_1.bin','r');
U = fread(f, 'float64');
fclose(f);
U = reshape(U,ns,nz,nx);


maxamp = max(max(max(U)));       % max amplitude calculation for clip
clip = 0;  % clip value
ampclip = (1-clip/100)*maxamp;
data    = U; 
data(data > ampclip)= ampclip;   % clipping positive
data(data < -ampclip)= -ampclip; % clipping negative


for is=1:ns
    imagesc(reshape(data(is,:,:),nz,nx));   
    colormap(redblue);
    %caxis([-6e-15 6e-15])
    caxis([-max(data(200,:,:),[],'all')  max(data(200,:,:),[],'all') ])

    colorbar

    pause(.0009) 
    
end

