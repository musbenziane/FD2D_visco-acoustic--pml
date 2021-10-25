clear; close all; clc;


disp('###### Reading Parameters ######')
fp = fopen("para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

%Parsing parameters
nt    = A{1}(1);
nx    = A{1}(2);
nz    = A{1}(3);
dx    = A{1}(4);
dz    = A{1}(5);
dt    = A{1}(6);
f0    = A{1}(7);
insrc = A{1}(8); % if 1 shot is slected, it will be at this position
isnap = A{1}(9);
nshot = A{1}(10);
shotz = A{1}(11);
nrec  = A{1}(12);
gWidth= A{1}(13);

nxs   = nx;
nx    = nx + 2 * gWidth;

sp    = floor(nxs/nrec);
rec   = gWidth:sp:nx;
sp    = floor(nxs/nshot);
shotp = gWidth:sp:nx - 2-gWidth;

% shot positions
if nshot==1
    isrc      = [shotz insrc]';
else
    isrc      = ones(2,nshot);
    isrc(1,:) = shotz;
    isrc(2,:) = shotp;
end

disp('##### Reading Model Files ######')
f     = fopen("c.bin","r");
c     = fread(f,"float32");
fclose(f);
c     = reshape(c,nz,nx);

cSize = [nz nx];
mSize = size(c);
if (mSize(1) ~= cSize(1) || mSize(2) ~= cSize(2))
    error_message = [ 'The model size does not match the computational domain size.\n' , ...
                  'Generate a new model.\n' , ...
                  'Program has been terminated' ];
    error( 'o:t' , error_message );
end

lam_min = min(min(c))/(f0*2.5);
gppw    = lam_min/dx;
tol     = 10;
U       = zeros(round(nt/isnap),nz,nx);
t       = 0:dt:nt*dt;
T       = 1/f0;

source  = (ricker(f0,nt,dt));

disp('########## CFL CHECK ##########');
cfl = (dt/dx)*max(max(c));
if (cfl>.9)
    error('Simulation unstable, decrease dt.');
else
    disp('Simulation is stable');
end

disp('# SPATIAL DISCRETEZATION CHECK #');

if (gppw<tol)
    error('Less than 8 grid points per wavelength are used, you risk numerical dispersion')
else
    disp('More than 8 grid points per wavelength are used, you are good');
end

u    = zeros(nz,nx);
uold = zeros(nz,nx);
dux  = u; duz = u;


% Absorbing Boundary Conditions.

g_z=ones(nz);
g_x=ones(nx);

  for k=1:gWidth
     att        = exp(-(0.06*(gWidth-k)/gWidth)^2);
     g_z(k)     = att;
     g_z(nz-k+1)= att;
     
     g_x(k)     = att;
     g_x(nx-k+1)= att;
  end
  
  gzx = ones(nz,nx);
  for i2=1:nx
     for i1=1:nz
        gzx(i1,i2) = g_z(i1) .* g_x(i2);
     end
  end

disp('####### Begin shot loop #######')

w1 = waitbar(0, 'Starting');
for is=1:length(shotp)
    

 
    w2 = waitbar(0, 'Starting');
    pos_w1=get(w1,'position');
    pos_w2=[pos_w1(1) pos_w1(2)+pos_w1(4) pos_w1(3) pos_w1(4)];
    set(w2,'position',pos_w2,'doublebuffer','on')
    waitbar(is/nshot, w1, sprintf('Shot loop progress: %d / %d', is,nshot));

    k = 1;
    shot_g =  zeros(nt,length(rec));  
    
    disp('####### Begin time loop #######')    
    for it=2:nt
        unew = zeros(nz,nx);


        unew(isrc(1,is),isrc(2,is)) = unew(isrc(1),isrc(2)) + dt^2 * source(it);

        for ix=2:nx-1     
            for iz=2:nz-1
                unew(iz,ix) = unew(iz,ix) + (c(iz,ix)^2) * (dt^2) *((((u(iz,ix+1) - 2 * u(iz,ix) + u(iz,ix-1)) / dz^2)) +...
                     ((u(iz+1,ix) - 2 * u(iz,ix) + u(iz-1,ix)) / dx^2)) + 2 * u(iz,ix) - uold(iz,ix); 
            end
        end
 
        uold = u    .* gzx;
        u    = unew .* gzx;
        
        
        for ir=1:length(rec)
            shot_g(it,ir) = unew(2,rec(ir));
        end
        if (mod(it,isnap)==0)

            U(k,:,:) = unew;
            k = k + 1;

        end
        waitbar(it/nt, w2, sprintf('Computation in progress: %d %%', floor(it/nt*100)));

    end
    
    filenameU  = strcat("OUTPUT/field_",mat2str(is),".bin");
    filenameSG = strcat("OUTPUT/seism_",mat2str(is),".bin");
    
    f1         = fopen(filenameU,"w");
    fwrite(f1,U,'float64');
    
    f2         = fopen(filenameSG,"w");
    fwrite(f2,shot_g,'float64');
    
    fclose(f1); fclose(f2);
    
end