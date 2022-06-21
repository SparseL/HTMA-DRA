function bestval = runcompe(fname, fun, D, NP, Max_Gen, runindex, fid);

if fun == 1 | fun==6 | fun==3 | fun==4
    XRmin = -100*ones(NP,D); 
    XRmax = 100*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 2
    XRmin = -10*ones(NP,D); 
    XRmax = 10*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 5
    XRmin = -30*ones(NP,D); 
    XRmax = 30*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 7
    XRmin = -1.28*ones(NP,D); 
    XRmax = 1.28*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 8
    XRmin = -500*ones(NP,D); 
    XRmax = 500*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 9
    XRmin = -5.12*ones(NP,D); 
    XRmax = 5.12*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 10
    XRmin = -32*ones(NP,D); 
    XRmax = 32*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 11
    XRmin = -600*ones(NP,D); 
    XRmax = 600*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun==12 | fun==13
    XRmin = -50*ones(NP,D); 
    XRmax = 50*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;   
end 
if fun == 14
    XRmin = -65.536*ones(NP,D); 
    XRmax = 65.536*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun==15 | fun==16
    XRmin = -5*ones(NP,D); 
    XRmax = 5*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;   
end
if fun == 17
    XRmin = ones(NP,1) * [-5  0];
    XRmax = ones(NP,1) * [10 15];
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun == 18
    XRmin = -2*ones(NP,D); 
    XRmax = 2*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;
end
if fun==19 | fun==20
    XRmin = 0*ones(NP,D); 
    XRmax = 1*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;   
end
if fun==21 | fun==22 | fun==23
    XRmin = 0*ones(NP,D); 
    XRmax = 10*ones(NP,D); 
    Lbound = XRmin;
    Ubound = XRmax;   
end

bestval = sansde(fname, fun, D, Lbound, Ubound, NP, Max_Gen, runindex, fid);

