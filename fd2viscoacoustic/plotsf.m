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
dx    = A{1}(4);
nx    = nx + npml*2;

fid      = fopen("Acqui/rcv.acqui","r");
format   = "%f %f";
sizeX    = [2 inf];
rcvpos     = fscanf(fid,format, sizeX);
fclose(fid);
ircv     = round((rcvpos(1,:))/dx+1); % rcv indices (1)
jrcv     = round((rcvpos(2,:))/dx+1); % rcv indices (2)
nrcv     = length(ircv);              % rcv number

f = fopen('OUTPUT/seis_1.bin','r');
U = fread(f, 'float64');
fclose(f);

U = reshape(U,nt,nrcv);


maxamp = max(max(max(U)));  % max amplitude calculation for clip
clip = 95;  % clip value
ampclip = (1-clip/100)*maxamp;
data    = U; 
data(data > ampclip)= ampclip; % clipping positive
data(data < -ampclip)= -ampclip; % clipping positive


imagesc(data)
colormap(redblue)
