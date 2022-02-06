clear; close all; clc;
try
    delete 'fd2dacoustic.log'
catch
end


diary 'fd2dacoustic.log'
diary on
disp('####################################################################')
disp('####                    Reading Parameters                      ####')
disp('####################################################################')

fp = fopen("para.in","r");
formatSpec = "%f %s";
A = textscan(fp,formatSpec);
fclose(fp); 

%Parsing parameters
nt       = A{1}(1);
nx       = A{1}(2);
nz       = A{1}(3);
dx       = A{1}(4);
dz       = A{1}(5);
dt       = A{1}(6);
f0       = A{1}(7);
insrc    = A{1}(8); % if 1 shot is slected, it will be at this position
isnap    = A{1}(9);
nshot    = A{1}(10);
shotz    = A{1}(11);
nrec     = A{1}(12);
irec     = A{1}(13);
zrec     = A{1}(14);
gWidth   = A{1}(15);
gWidthz  = A{1}(16);
attConst = A{1}(17);

sp    = floor(nx/nshot);
shotp = gWidth:sp:nx-gWidth;

% shot positions
if nshot==1
    isrc      = [shotz insrc]';
else
    isrc      = ones(2,nshot);
    isrc(1,:) = shotz;
    isrc(2,:) = shotp;
end

% rec positions
if nrec==1
    rec   = irec;
else
    sp    = floor((nx-2*gWidth)/(nrec));
    rec   = gWidth+1:sp:nx-gWidth-1;
end

fprintf('\n');
disp('####################################################################')
disp('####                  Simulation parameters:                    ####')
fprintf('####                  nt            ->   %d                          \n',nt);
fprintf('####                  dt            ->   %1.2E                       \n',dt);
fprintf('####                  nx            ->   %d                          \n',nx);
fprintf('####                  nz            ->   %d                          \n',nz);
fprintf('####                  nx            ->   %d                          \n',nx);
fprintf('####                  dz            ->   %d                          \n',nz);
fprintf('####                  dx            ->   %d                          \n',nx);
fprintf('####                  Ricker f0     ->   %d                          \n',f0);
fprintf('####         &        Sponge size   ->   %d                          \n',gWidth);
fprintf('####                  Shots num     ->   %d                          \n',nshot);
fprintf('####                  Receiver num  ->   %d                          \n',nrec);
fprintf('####                  Snap interval ->   %d                          \n',isnap);
disp('####################################################################')
fprintf('\n');

disp('####################################################################')
disp('####                    Reading Model Files                     ####')
disp('####################################################################')

f     = fopen("ctrial.bin","r");
c     = fread(f,"float32");
fclose(f);


cSize = [nz nx];
mSize = size(c);

try
    c     = reshape(c,nz,nx);
catch
    error_message = [ 'The model size does not match the computational domain size.\n' , ...
                  'Generate a new model.\n' , ...
                  'Program has been terminated' ];
    error( 'o:t' , error_message );

end
lam_min = min(min(c))/(f0*2.5);
gppw    = lam_min/dx;
tol     = 9;
U       = zeros(round(nt/isnap),nz,nx);
t       = 0:dt:nt*dt;
T       = 1/f0;

source  = (ricker(f0,nt,dt));


fprintf('\n');
disp('####################################################################')
disp('####                     CFL Creterion check                    ####')

cfl = (dt/dx)*max(max(c));
if (cfl>.9)
    error('####             Simulation unstable, decrease dt.         ####');
else
    disp('####                     Simulation is stable               ####')
 fprintf('####                     Courant number ->   %1.2f          #### \n',cfl);

end
disp('####################################################################')

fprintf('\n');
disp('####################################################################')
disp('####                SPATIAL DISCRETEZATION CHECK                ####')
if (gppw<tol)
    error('####    Numerical dispersion might be present, decrease dx|dz   ####')
else
    disp('####                Spatial sampling seems fine                 ####')
end
disp('####################################################################')


% Absorbing Boundary Conditions.

g_z      = ones(nz);
g_x      = ones(nx);

  for k=1:gWidth
     att        = exp(-(attConst*(gWidth-k)/gWidth)^2);
     g_z(nz-k+1)= att;
     g_x(k)     = att;
     g_x(nx-k+1)= att;
  end
  
  for k=1:gWidthz
     att        = exp(-(attConst*(gWidthz-k)/20)^2);
     g_z(k)     = att;
  end
  
  gzx = ones(nz,nx);
  for ix=1:nx
     for iz=1:nz
        gzx(iz,ix) = g_z(iz) .* g_x(ix);
     end
  end
fprintf('\n');
disp('####################################################################')
disp('####                      Begin shot loop                       ####')
disp('####################################################################')

w1 = waitbar(0, 'Starting');

    fprintf('\n');


for is=1:nshot

 
    disp('####################################################################')
    fprintf('####                      Shot number %d                        ####\n',is);
    disp('####################################################################')
    w2 = waitbar(0, 'Starting');
    pos_w1=get(w1,'position');
    pos_w2=[pos_w1(1) pos_w1(2)+pos_w1(4) pos_w1(3) pos_w1(4)];
    set(w2,'position',pos_w2,'doublebuffer','on')
    waitbar(is/nshot, w1, sprintf('Shot loop progress: %d / %d', is,nshot));

    k = 1;
    shot_g =  zeros(nt,nrec);  
    u    = zeros(nz,nx);
    uold = zeros(nz,nx);
    dux  = u; duz = u;

    %####################################################################
    %####                      Begin time loop                       ####
    %####################################################################
  
    
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
        
        
        for ir=1:nrec
            shot_g(it,ir) = unew(zrec,rec(ir));
        end
        if (mod(it,isnap)==0)

            U(k,:,:) = unew;
            k = k + 1;

        end
        waitbar(it/nt, w2, sprintf('Computation in progress: %d %%', floor(it/nt*100)));

    end
    
    
    %####################################################################
    %####                      Output snapshot                       ####
    %####                      & seismograms                         ####
    %####################################################################
    
    filenameU  = strcat("OUTPUT/field_",mat2str(is),".bin");
    filenameSG = strcat("OUTPUT/seis_",mat2str(is),".bin");
    
    f1         = fopen(filenameU,"w");
    fwrite(f1,U,'float64');
    
    f2         = fopen(filenameSG,"w");
    fwrite(f2,shot_g,'float64');
    
    fclose(f1); fclose(f2);
    
end

fprintf('\n');

disp('####################################################################')
disp('####                     Program is ending                     ####')
disp('####                   files are in OUTPUT/                    ####')
disp('####                  log file: fd2dacoustic.log               ####')
disp('####################################################################')


diary off
