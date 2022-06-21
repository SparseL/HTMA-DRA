
clear;
NP=100;
runs=50;
global initial_flag;

dim = ones(1,13)*30;
dim = [dim 2 4 2 2 2 4 6 4 4 4];

Gen = ones(1, 13)*1500;
Gen(5) = 5000;
Gen = [Gen ones(1,10)*200];
Gen(15) = 1500;

myfunc = 1:23;
addpath('benchmark');
for fun=1:23
    func_num = myfunc(fun);
    D = dim(func_num);
    Max_Gen = Gen(func_num);
    bestval = [];
    for runindex=1:runs
        filename = sprintf('trace/tracef%02d_%02d.txt', func_num, runindex);
        fid = fopen(filename, 'w');
        initial_flag = 0;
        val = runcompe('benchmark_func', func_num, D, NP, Max_Gen, runindex, fid);
        bestval = [bestval; val];
        fclose(fid);
    end
    
    filename = sprintf('result/bestf%02d.txt', func_num);
    fid = fopen(filename, 'w');
    fprintf(fid, '%g\n', bestval);
    fclose(fid);
    
    filename = sprintf('result/meanf%02d.txt', func_num);
    fid = fopen(filename, 'w');
    fprintf(fid, '%g\n', mean(bestval));
    fclose(fid);
    
    filename = sprintf('result/stdf%02d.txt', func_num);
    fid = fopen(filename, 'w');
    fprintf(fid, '%g\n', std(bestval));
    fclose(fid);
end

