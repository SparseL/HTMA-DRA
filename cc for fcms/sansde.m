
function bestval = sansde(fname, fun, D, Lbound, Ubound, NP, itermax, runindex, fid)

refresh = 500;
rand('state',sum(100*clock));

CRstart = 0.5;

tracerst    = [];
traceccm    = [];
tracelinkp  = [];
traceflinkp = [];

ccm = CRstart;

%initialize pop
pop = zeros(NP,D);

F = zeros(NP,1);

pop = Lbound + rand(NP, D) .* (Ubound-Lbound);
if fun == 7
    XRmin = 0*ones(NP,D); 
    XRmax = 600*ones(NP,D);
    pop = XRmin + rand(NP, D) .* (XRmax-XRmin);
    XRmin = [];
    XRmax = [];
end
if fun == 25
    XRmin = -2*ones(NP,D); 
    XRmax = 5*ones(NP,D); 
    pop = XRmin + rand(NP, D) .* (XRmax-XRmin);
    XRmin = [];
    XRmax = [];
end

popold    = zeros(size(pop));     % toggle population
val       = zeros(NP, 1);         

linkp = 0.5;
l1 = 1; l2 = 1; nl1 = 1; nl2 = 1;

fp = 0.5;
ns1 = 1; nf1 = 1; ns2 = 1; nf2 = 1;

%------Evaluate the best member after initialization----------------------
val = feval(fname, pop, fun);
[bestval ibest] = min(val);
bestmem = pop(ibest, :);

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize DE_gbestber  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);

iter = 0;
while iter < itermax
    popold = pop;                   % save the old population
    
    ind = randperm(4);              % index pointer array
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);                
    rt = rem(rot+ind(3),NP);
    a4  = a3(rt+1);               
    rt = rem(rot+ind(4),NP);
    a5  = a4(rt+1); 
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    bm = ones(NP, 1) * bestmem;
    
    if rem(iter,25)==0
        if (iter~=0) & (~isempty(cc_rec))
           ccm = sum(f_rec.*cc_rec)/sum(f_rec);
        end
        cc_rec = [];
        f_rec = [];
    end
    traceccm = [traceccm; ccm];
    
    if rem(iter,5)==0
        cc = normrnd(ccm, 0.1, NP*3, 1);
        index = find((cc < 1) & (cc > 0));
        cc = cc(index(1:NP));
    end

    fst1 = (rand(NP,1) <= fp);
    fst2 = 1-fst1;
    traceflinkp = [traceflinkp; fp sum(fst1) sum(fst2)];
    
    fst1_index = find(fst1 ~= 0);
    fst2_index = find(fst1 == 0);

    tmp = normrnd(0.5, 0.3, NP, 1);
    F(fst1_index) = tmp(fst1_index);
    
    tmp = normrnd(0, 1, NP, 1) ./ normrnd(0, 1, NP, 1);
    F(fst2_index) = tmp(fst2_index);
    
    F = abs(F);
    
    %all random numbers < CR are 1, 0 otherwise
    aa = rand(NP,D) < repmat(cc,1,D); 
    index = find(sum(aa') == 0);
    tmpsize = size(index, 2);
    for k=1:tmpsize
        bb = ceil(D*rand);
        aa(index(k), bb) = 1;
    end
        
    mui = aa;
    mpo = mui < 0.5;                % inverse mask to mui
    
    st1 = (rand(NP,1) <= linkp);
    st2 = 1-st1;
    tracelinkp = [tracelinkp; linkp sum(st1) sum(st2)];
    
    st1_index = find(st1 ~= 0);
    st2_index = find(st1 == 0);
    
    if ~isempty(st1_index)
        ui(st1_index,:) = pm1(st1_index,:) + repmat(F(st1_index,:),1,D).*(pm2(st1_index,:) - pm3(st1_index,:));
    
        ui(st1_index,:) = popold(st1_index,:).*mpo(st1_index,:) + ui(st1_index,:).*mui(st1_index,:);
    end
   
    if ~isempty(st2_index)
        ui(st2_index,:) = popold(st2_index,:)+repmat(F(st2_index,:),1,D).*(bm(st2_index,:)-popold(st2_index,:)) + repmat(F(st2_index,:),1,D).*(pm1(st2_index,:) - pm2(st2_index,:) + pm3(st2_index,:) - pm4(st2_index,:));
        
        ui(st2_index,:) = popold(st2_index,:).*mpo(st2_index,:) + ui(st2_index,:).*mui(st2_index,:);
    end
    

    index = find(ui > Ubound);
    ui(index) = Ubound(index) - mod((ui(index)-Ubound(index)), (Ubound(index)-Lbound(index)));
    index = find(ui < Lbound);
    ui(index) = Lbound(index) + mod((Lbound(index)-ui(index)), (Ubound(index)-Lbound(index)));

    tempval = feval(fname, ui, fun);
    for i=1:NP
       if (tempval(i) <= val(i))  
            if (tempval(i) < val(i))
               cc_rec = [cc_rec cc(i,1)];
               f_rec = [f_rec (val(i) - tempval(i))];
            end
            
            pop(i,:) = ui(i,:);  
            val(i)   = tempval(i);  
            
            l1 = l1 + st1(i);
            l2 = l2 + st2(i);

            ns1 = ns1 + fst1(i);
            ns2 = ns2 + fst2(i);

        else
            nl1 = nl1 + st1(i);
            nl2 = nl2 + st2(i);

            nf1 = nf1 + fst1(i);
            nf2 = nf2 + fst2(i);
        end
    end
 
    if (rem(iter,50) == 0) & (iter~=0);
        linkp = (l1/(l1+nl1))/(l1/(l1+nl1)+l2/(l2+nl2));
        l1 = 1; l2 = 1; nl1 = 1; nl2 = 1;
        
        fp = (ns1 * (ns2 + nf2))/(ns2 * (ns1 + nf1) + ns1 * (ns2 + nf2));
        ns1 = 1; nf1 = 1; ns2 = 1; nf2 = 1;
    end

    [bestval, ibest] = min(val);
    bestmem = pop(ibest, :);
    tracerst = [tracerst; bestval];
       
    iter = iter + 1;
    
    if (rem(iter, 1000) == 0)
      fprintf(fid, '%g\n', tracerst);
      tracerst = [];
  end
                             
  if (refresh > 0)
    if (rem(iter, refresh) == 0)
        fprintf(1, 'linkp = %f, fp = %f, ccm = %g\n', linkp, fp, ccm);
        fprintf(1, 'fun = %d, run = %d, Iter = %d, Dim = %d, NP = %d, bestval = %g\n\n', fun, runindex, iter, D, NP, bestval);
    end
  end
end 

 if (~isempty(tracerst))
    fprintf(fid, '%g\n', tracerst);
end

fname = sprintf('runtime/ccm_f%02d_%02d.txt', fun, runindex);
ccmpid = fopen(fname, 'w');
fprintf(ccmpid, '%g\n', traceccm);
fclose(ccmpid);

fname = sprintf('runtime/linkp_f%02d_%02d.txt', fun, runindex);
linkpid = fopen(fname, 'w');
for i = 1:size(tracelinkp)
   fprintf(linkpid, '%g ',  tracelinkp(i, 1));
   fprintf(linkpid, '%d ',  tracelinkp(i, 2));
   fprintf(linkpid, '%d\n', tracelinkp(i, 3));
end
fclose(linkpid);

fname = sprintf('runtime/flinkp_f%02d_%02d.txt', fun, runindex);
flinkpid = fopen(fname, 'w');
for i = 1:size(traceflinkp)
   fprintf(flinkpid, '%g ',  traceflinkp(i, 1));
   fprintf(flinkpid, '%d ',  traceflinkp(i, 2));
   fprintf(flinkpid, '%d\n', traceflinkp(i, 3));
end
fclose(flinkpid);

