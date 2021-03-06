function source = ricker(f0,nt,dt)    

        source   = 0;
        temp     = 0;
        pt       = 1 / f0;
        t0       = pt/dt;
        a_ricker = 4./pt;

        for it=1:nt+1
            t=(it-t0)*dt;
            temp(it) = -2.*a_ricker*t*exp(-(a_ricker*t)^2);
        end
        
        
        for it=1:nt
            source(it) = temp(it+1) - temp(it);
            source(it) = -source(it);
        end 

end