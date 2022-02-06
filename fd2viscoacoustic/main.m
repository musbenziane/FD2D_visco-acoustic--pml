clear; close all; clc;

try
    delete 'fd2dviscoacoustic.log'
catch
end


diary 'fd2dviscoacoustic.log'
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
dt       = A{1}(5);
f0       = A{1}(6);
npml     = A{1}(7);
is_pml   = A{1}(8);
is_FS    = A{1}(8);
isnap    = A{1}(10);


%####################################################################
%####                Init. src, time, position                   ####
%####################################################################
t       = 0:dt:(nt-1)*dt;
source  = (fn_ricker(f0,nt,dt));
x       = 0:dx:(nx-1)*dx;
z       = 0:dx:(nz-1)*dx;

%####################################################################
%####                Read acquisition geometry                   ####
%####################################################################
fid      = fopen("Acqui/rcv.acqui","r");
format   = "%f %f";
sizeX    = [2 inf];
rcvpos     = fscanf(fid,format, sizeX);
fclose(fid);

fid      = fopen("Acqui/src.acqui","r");
format   = "%f %f";
sizeX    = [2 inf];
srcpos   = fscanf(fid,format, sizeX);
fclose(fid);

isrc     = round((srcpos(1,:))/dx+1); % src indices (1)
jsrc     = round((srcpos(2,:))/dx+1); % src indices (2)
nshot    = length(isrc);              % src number

ircv     = round((rcvpos(1,:))/dx+1); % rcv indices (1)
jrcv     = round((rcvpos(2,:))/dx+1); % rcv indices (2)
nrcv     = length(ircv);              % rcv number


fprintf('\n');
disp('####################################################################')
disp('####                  Simulation parameters:                    ####')
fprintf('####                  nt            ->   %d                          \n',nt);
fprintf('####                  dt            ->   %1.2E                       \n',dt);
fprintf('####                  nx            ->   %d                          \n',nx);
fprintf('####                  nz            ->   %d                          \n',nz);
fprintf('####                  nx            ->   %d                          \n',nx);
fprintf('####                  Ricker f0     ->   %d                          \n',f0);
fprintf('####         &        PML points    ->   %d                          \n',npml);
fprintf('####                  Shots num     ->   %d                          \n',nshot);
fprintf('####                  Receiver num  ->   %d                          \n',nrcv);
fprintf('####                  Snap interval ->   %d                          \n',isnap);
disp('####################################################################')
fprintf('\n');

disp('####################################################################')
disp('####                    Reading Model Files                     ####')
disp('####################################################################')

f         = fopen("Mod_files/c.bin","r");
c         = fread(f,"float32");
fclose(f);

f         = fopen("Mod_files/rho.bin","r");
rho       = fread(f,"float32");
fclose(f);

f         = fopen("Mod_files/qf.bin","r");
qf        = fread(f,"float32");
fclose(f);


try
    c     = reshape(c,nz,nx);
    rho   = reshape(rho,nz,nx);
    qf    = reshape(qf,nz,nx);

catch
    error_message = [ 'The model size does not match the computational domain size.\n' , ...
                  'Generate a new model.\n' , ...
                  'Program has been terminated' ];
    error( 'o:t' , error_message );

end

%####              Extend model for PML                ####
if is_pml
    c   = [repmat(c(:,nz),1,npml) c repmat(c(:,1),1,npml)]; 
    rho = [repmat(rho(:,nz),1,npml) rho repmat(rho(:,1),1,npml)];
    qf  = [repmat(qf(:,nz),1,npml) qf repmat(qf(:,1),1,npml)];

    c   = [c; repmat(c(nz,:),npml,1)];
    rho = [rho; repmat(rho(nz,:),npml,1)];
    qf  = [qf; repmat(qf(nz,:),npml,1)];
    nx = nx + 2*npml;
    nz = nz + npml;
end

K       = c .^2 .* rho;                           % Relaxed elastic modulus
tau_sigma = (sqrt(1+1.0./ qf.^2)-1.0./qf)./f0;    % Relaxation time: stress
tau_epsi  = 1.0./(f0^2.0.*tau_sigma);             % Relaxation time: strain


fprintf('\n');
disp('####################################################################')
disp('####                     CFL Creterion check                    ####')

cfl = (dt/dx)*max(max(c));
if (cfl>.34)
    error('####             Simulation unstable, decrease dt.         ####');
else
    disp('####                     Simulation is stable               ####')
 fprintf('####                     Courant number ->   %1.2f          #### \n',cfl);

end
disp('####################################################################')

fprintf('\n');
disp('####################################################################')
disp('####                SPATIAL DISCRETEZATION CHECK                ####')
lam_min = min(min(c))/(f0*2.5);
gppw    = lam_min/dx;
tol     = 12;
if (gppw<tol)
    error('####    Numerical dispersion might be present, decrease dx|dz   ####')
else
    disp('####                Spatial sampling seems fine                 ####')
end
disp('####################################################################')

                                       

%####################################################################
%####                   Stencils' coeffecient                    ####
%####################################################################
%{   
@article{article,
author = {Liu, Yang and Sen, Mrinal},
year = {2009},
month = {10},
pages = {459-474},
title = {An implicit staggered-grid finite-difference method for seismic modeling},
volume = {179},
journal = {Geophysical Journal International - GEOPHYS J INT},
doi = {10.1111/j.1365-246X.2009.04305.x}
}
%}
% Program is set up for 8th, can be switched to 4th easily. 
c1=1225.0/1024.0; 
c2=-245.0/3072.0;
c3=49.0/5120.0;   
c4=-5.0/7168.0;    
%c1=9.0/8.0; c2=-1.0/24.0;             % 4th 


%####################################################################
%####                           PML    Profile                   ####
%####################################################################

sigma_m = 450;
wPML    = dx*(npml-1);
tmp     = zeros(npml,1);
sigmaZ  = zeros(nz,1);
sigmaX  = zeros(nx,1);

for ix=1:npml
  sm          = real(ix-1)*dx/wPML;
  tmp(ix) = sigma_m*sm.^2;
end

sigmaZ(nz-npml+1:end) = tmp;
sigmaX(nx-npml+1:end) = tmp;
sigmaX(1:npml)        = flipud(tmp);

if not(is_pml)
    sigmaX(:) = 0;
    sigmaZ(:) = 0;
end


fprintf('\n');
disp('####################################################################')
disp('####                      Begin shot loop                       ####')
disp('####################################################################')

fprintf('\n');

for is=1:nshot

 
    disp('####################################################################')
    fprintf('####                      Shot number %d                        ####\n',is);
    disp('####################################################################')
    
    w1 = waitbar(0, 'Starting');

    k      = 1;
    seis_p =  zeros(nt,nrcv);  
    seis_u = seis_p;
    seis_w = seis_p;
    u      = zeros(nz,nx);
    p      = u; 
    w      = u;
    r      = u;
    px     = p;
    pz     = p;
    U      = zeros(round(nt/isnap),nz,nx);
    W      = zeros(round(nt/isnap),nz,nx);
    P      = zeros(round(nt/isnap),nz,nx);

    %####################################################################
    %####                      Begin time loop                       ####
    %####################################################################
  
    
    for it=2:nt

        
        for iz=5:nz-4
            for ix=5:nx-4

                temp      = K(iz,ix)*dt*(tau_epsi(iz,ix)/tau_sigma(iz,ix)-1.0)/dx/tau_sigma(iz,ix);
                r(iz,ix)  =  (1.0-dt/tau_sigma(iz,ix))*r(iz,ix)-temp          * ...
                            (c1*(u(iz,ix+1)-u(iz,ix)+w(iz,ix)-w(iz-1,ix))     + ...
                             c2*(u(iz,ix+2)-u(iz,ix-1)+w(iz+1,ix)-w(iz-2,ix)) + ...
                             c3*(u(iz,ix+3)-u(iz,ix-2)+w(iz+2,ix)-w(iz-3,ix)) + ...
                             c4*(u(iz,ix+4)-u(iz,ix-3)+w(iz+3,ix)-w(iz-4,ix)));
                      
                temp      = K(iz,ix)*dt*tau_epsi(iz,ix)/dx/tau_sigma(iz,ix);
                pz(iz,ix) = (1-dt*sigmaZ(iz)) * pz(iz,ix)-r(iz,ix)*dt-temp          * ...    
                            (c1*(w(iz,ix)-w(iz-1,ix))  + c2*(w(iz+1,ix)-w(iz-2,ix)) + ...
                             c3*(w(iz+2,ix)-w(iz-3,ix))+ c4*(w(iz+3,ix)-w(iz-4,ix)));
            
                px(iz,ix) = (1-dt*sigmaX(ix)) * px(iz,ix)-r(iz,ix)*dt-temp         * ...    
                            (c1*(u(iz,ix+1)-u(iz,ix))  + c2*(u(iz,ix+2)-u(iz,ix-1))+ ...
                             c3*(u(iz,ix+3)-u(iz,ix-2))+ c4*(u(iz,ix+4)-u(iz,ix-3)));

                p(iz,ix) = pz(iz,ix) + px(iz,ix); 
                  
                temp     = dt/(rho(iz,ix)*dx);
                w(iz,ix) = (1-dt*sigmaZ(iz)) * w(iz,ix)-temp   * ...
                                    (c1*(p(iz+1,ix)-p(iz,ix))  +...
                                     c2*(p(iz+2,ix)-p(iz-1,ix))+...  
                                     c3*(p(iz+3,ix)-p(iz-2,ix))+...
                                     c4*(p(iz+4,ix)-p(iz-3,ix)));
                  
                  
                temp     = dt/(rho(iz,ix)*dx);
                u(iz,ix) = (1-dt*sigmaX(ix)) * u(iz,ix) - temp  * ...
                                    (c1*(p(iz,ix)-p(iz,ix-1))   +...
                                     c2*(p(iz,ix+1)-p(iz,ix-2)) +...  
                                     c3*(p(iz,ix+2)-p(iz,ix-3)) +...
                                     c4*(p(iz,ix+3)-p(iz,ix-4)));

            end
        end
        
        p(isrc(is),jsrc(is)) = p(isrc(is),isrc(is)) + dt * source(it);
        
        if is_FS
            p(5,:) = 0;
            pz(5,:)= 0;
            px(5,:)= 0;
        end
        
        for ir=1:nrcv
            seis_p(it,ir) = p(ircv(ir),jrcv(ir));
            seis_w(it,ir) = w(ircv(ir),jrcv(ir));
            seis_u(it,ir) = u(ircv(ir),jrcv(ir));
        end
        
        if (mod(it,isnap)==0)
            P(k,:,:) = p;
            W(k,:,:) = w;
            U(k,:,:) = u;
            k = k + 1;
        end
        waitbar(it/nt, w1, sprintf('Computation: %d %% | shot %d of %d ', floor(it/nt*100),is,nshot));
    end
    
    
    %####################################################################
    %####                      Output snapshot                       ####
    %####                      & seismograms                         ####
    %####################################################################
   
    filenameU  = strcat("OUTPUT/field_P_",mat2str(is),".bin");
    filenameSG = strcat("OUTPUT/seis_P_",mat2str(is),".bin");
    
    f1         = fopen(filenameU,"w");
    fwrite(f1,P,'float64');
    
    f2         = fopen(filenameSG,"w");
    fwrite(f2,seis_p,'float64');
    
    fclose(f1); fclose(f2);
    
    
    filenameU  = strcat("OUTPUT/field_W_",mat2str(is),".bin");
    filenameSG = strcat("OUTPUT/seis_W_",mat2str(is),".bin");
    
    f1         = fopen(filenameU,"w");
    fwrite(f1,W,'float64');
    
    f2         = fopen(filenameSG,"w");
    fwrite(f2,seis_w,'float64');
    
    fclose(f1); fclose(f2);
    
    filenameU  = strcat("OUTPUT/field_U_",mat2str(is),".bin");
    filenameSG = strcat("OUTPUT/seis_U_",mat2str(is),".bin");
    
    f1         = fopen(filenameU,"w");
    fwrite(f1,U,'float64');
    
    f2         = fopen(filenameSG,"w");
    fwrite(f2,seis_u,'float64');
    
    fclose(f1); fclose(f2);
    
end

fprintf('\n');

disp('####################################################################')
disp('####                     Program is ending                      ####')
disp('####                   files are in OUTPUT/                     ####')
disp('####                  log file: fd2dacoustic.log                ####')
disp('####################################################################')


diary off