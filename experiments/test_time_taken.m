%% Set up problem

addpath('code')
addpath('utils')

m = 5e6;
sizes = round(logspace(1,3,5));

trials = 50; 

cond_A = 10^4;
res_size = 10^-3;

theta = Inf;

qrs = zeros(length(sizes), 1);
wedins = zeros(length(sizes), 1);
itsks = zeros(length(sizes), 2);
sps = zeros(length(sizes), 2);

real_run = true;

for i = 1:5
    A = [];
    [A,b,x,r] = create_ls_problem(m,sizes(i),cond_A,res_size);

    % QR
    tic; x_qr = A\b; qrs(i) = toc;
    qrs(i)

    % Iterative sketching
    tic; x_is = iterative_sketching(A,b,20*sizes(i)); itsks(i,1) = toc;
    itsks(i,2) = norm(x_is - x_qr) / norm(x_qr);
    itsks(i,1)

    % Sketch and precondition
    tic; x_sp = sketch_and_precondition(A,b,20*sizes(i),100); sps(i,1) = toc;
    sps(i,2) = norm(x_sp - x_qr) / norm(x_qr);
    sps(i,1)

    % Wedin
    [~,R] = qr(A,'econ');
    r_qr = b - A*x_qr;
    cond_A = cond(R)
    norm_A = norm(R)
    wedins(i) = 2.23 * cond_A * (1 + cond_A / norm_A / norm(x_qr) * norm(r_qr)) * eps;
end

if real_run
    save('data/results_timing.mat', 'qrs', 'itsks', 'sps', 'wedins')
end

%% Plot

close all
figure(1)
loglog(sizes, qrs(:,1), 'k', 'LineWidth', 4); hold on
loglog(sizes, itsks(:,1), ':', 'LineWidth',4,'Color',"#EDB120");
loglog(sizes, sps(:,1), ':', 'LineWidth',4,'Color',"#0072BD");
xlabel('Number of columns $n$')
ylabel('Time (sec)')
legend({'\texttt{mldivide}', 'Iterative sketching', 'Sketch and precondition'},'Location','best')

if real_run
    saveas(gcf, 'figs/sketch_times_taken.fig')
    saveas(gcf, 'figs/sketch_times_taken.png')
end

figure(2)
loglog(sizes, wedins, ':', 'LineWidth', 3, 'Color', "#A2142F"); hold on
loglog(sizes, itsks(:,2), ':', 'LineWidth',4,'Color',"#EDB120");
loglog(sizes, sps(:,2), ':', 'LineWidth',4,'Color',"#0072BD");
xlabel('Number of columns $n$')
ylabel('Forward error $\|\mbox{\boldmath $\widehat{x}$}_{\mbox{dir}}-\mbox{\boldmath $\widehat{x}$}\|/\|\mbox{\boldmath $\widehat{x}$}_{\mbox{dir}}\|$')

if real_run
    saveas(gcf, 'figs/sketch_errors.fig')
    saveas(gcf, 'figs/sketch_errors.png')
end