clear;
clc;


% WARNING all input files must be of the type reps=1


%     'PGDM'
%     'FISTA'
%     'SPGD'
%     'DIA'
%     'CVX'
% D = [2,4,8,16,32,64,128,256];
D = [2,4,8,16,32,64,128,256];

% TODO accept any number of repitions, and collate
PGDM = [nan];
FISTA = [nan];
SPGD = [nan];
DIA = [nan];
CVX = [nan];
for n = 1:6
    d = 2^n;
    DIR = ['gk_results_double_precision_single_thread_tinis/d',num2str(d),'/'];
    dirData = dir(DIR);
    dirData = dirData(~[dirData.isdir]);  %# Use only the file data
    fileNames = {dirData.name};

    disp(length(fileNames))
    for i = 1:length(fileNames)
        
        cd(DIR)
        file = fileNames(i);

        load(file{:,1})
        cd ../..
        l =  length(TIMETAKEN(1,:));
        PGDM(n,end+1:end+l) = TIMETAKEN(1,:);
        FISTA(n,end+1:end+l) = TIMETAKEN(2,:);
        SPGD(n,end+1:end+l) = TIMETAKEN(3,:);
        if length(TIMETAKEN(:,1))>3
            DIA(n,end+1:end+l) = TIMETAKEN(4,:);
            CVX(n,end+1:end+l) = TIMETAKEN(5,:);
        end
    end
end

PGDMvar = nanvar(PGDM')';
PGDMmean = nanmean(PGDM,2);
FISTAvar = nanvar(FISTA')';
FISTAmean = nanmean(FISTA,2);
SPGDvar = nanvar(SPGD')';
SPGDmean = nanmean(SPGD,2);
DIAvar = nanvar(DIA')';
DIAmean = nanmean(DIA,2);
CVXvar = nanvar(CVX')';
CVXmean = nanmean(CVX,2);
        
nnz(PGDM(1,:))
nnz(PGDM(2,:))
% t{1} = load('Result_d=2_reps=50');
% t{2} = load('Result_d=4_reps=40');
% t{3} = load('Result_d=8_reps=30');
% t{4} = load('Result_d=16_reps=20');
% t{5} = load('Result_d=32_reps=10');
% t{6} = load('Result_d=64_reps=4');
% t{7} = load('Result_d=128_reps=2');

% for i=1:5
%     PGDM(i) = mean(t{i}.TIMETAKEN(1,:));
%      FISTA(i) = mean(t{i}.TIMETAKEN(2,:));
%       SPGD(i) = mean(t{i}.TIMETAKEN(3,:));
%        DIA(i) = mean(t{i}.TIMETAKEN(4,:));
%         CVX(i) = mean(t{i}.TIMETAKEN(5,:));
% end
% for i=6:7
%     PGDM(i) = mean(t{i}.TIMETAKEN(1,:));
%      FISTA(i) = mean(t{i}.TIMETAKEN(2,:));
%       SPGD(i) = mean(t{i}.TIMETAKEN(3,:));
% end

figure(1)
hold off
x=[1 2 3 4 5 6 7 8];
x=[1 2 3 4 5 6];
errorbar(x,PGDMmean,PGDMvar);
hold on
errorbar(x,FISTAmean,FISTAvar)
errorbar(x,SPGDmean,SPGDvar)
errorbar([1 2 3 4 5 6 ],DIAmean,DIAvar)
errorbar([1 2 3 4 5 6 ],CVXmean,CVXvar)
title('double precision data from tinis')
set(gca,'yscale','log')
% semilogy([1 2 3 4 5 ],DIA)
% semilogy([1 2 3 4 5 ],CVX)
xlabel('Number of qubits')
ylabel('Time (s)')
% ylim([0.0001,1000])
ylim auto
legend('PGDM','FISTA','PGDB','DIA','CVX','Location','southeast')
saveas(gcf,'gk_DOUBLE_TINIS.png')