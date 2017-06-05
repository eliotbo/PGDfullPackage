%  function makeTimeGraph()

clc;clear;

if ~ismac && ~ispc
    cd cvx_linux
    cvx_setup
    cd ..
end

%add subfolders as paths
currentFolderContents = dir(pwd);
currentFolderContents (~[currentFolderContents.isdir]) = [];
for i = 3:length(currentFolderContents)
    addpath(['./' currentFolderContents(i).name]) ;
end

D=[4];
REPEAT = [ 1];
k=0;
for d=D
    rng('shuffle')
    param.d = d;
    param.purity = 0.5;             %purity of the density matrix
    param.counts = param.d*1E4;     %multication factor r
    param.theta = 90*pi/180;       %only used if d=2^n (many qubits).

    [A] = makeQubitILLMUB(log2(param.d ),param.theta);
    whos('A')
    if ispc
        memory
    end
    k=k+1;
    repeat = REPEAT(k);
    a=0;
    for reps = 1:repeat
        a=a+1;
        [fid,timeTaken] = example1(param,A);
        
        fields = fieldnames(timeTaken);
        for i = 1:numel(fields)
            FID(i,a) = abs(fid.(fields{i}));
            TIMETAKEN(i,a) = timeTaken.(fields{i});
        end
    end
    
    disp(fields)
    meanFid = 1-mean(FID,2)
    meanTime = mean(TIMETAKEN,2)
    stdFid = std(FID')'
    stdTime = std(TIMETAKEN')'
    
    save(['ResultT_d=' num2str(d) '_reps=' num2str(repeat) '_rand' num2str(round(rand(1)*10000000))],'FID','TIMETAKEN','fields')
    
end