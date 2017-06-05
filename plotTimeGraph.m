clear;
clc;

%     'PGDM'
%     'FISTA'
%     'SPGD'
%     'DIA'
%     'CVX'
    
t{1} = load('Result_d=2_reps=50');
t{2} = load('Result_d=4_reps=40');
t{3} = load('Result_d=8_reps=30');
t{4} = load('Result_d=16_reps=20');
t{5} = load('Result_d=32_reps=10');
t{6} = load('Result_d=64_reps=4');
t{7} = load('Result_d=128_reps=2');

for i=1:5
    PGDM(i) = mean(t{i}.TIMETAKEN(1,:));
     FISTA(i) = mean(t{i}.TIMETAKEN(2,:));
      SPGD(i) = mean(t{i}.TIMETAKEN(3,:));
       DIA(i) = mean(t{i}.TIMETAKEN(4,:));
        CVX(i) = mean(t{i}.TIMETAKEN(5,:));
end
for i=6:7
    PGDM(i) = mean(t{i}.TIMETAKEN(1,:));
     FISTA(i) = mean(t{i}.TIMETAKEN(2,:));
      SPGD(i) = mean(t{i}.TIMETAKEN(3,:));
end

figure(1)
hold off
x=[1 2 3 4 5 6 7];
semilogy(x,PGDM);
hold on
semilogy(x,FISTA)
semilogy(x,SPGD)
semilogy([1 2 3 4 5 ],DIA)
semilogy([1 2 3 4 5 ],CVX)
xlabel('Number of qubits')
ylabel('Time (s)')