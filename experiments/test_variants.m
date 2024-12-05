%% Set up problem

addpath('code')
addpath('utils')

m = 2000;
n = 50;

trials = 50; 

cond_A = 10^10;
res_size = 10^-6;

[A,b,x,r] = create_ls_problem(m,n,cond_A,res_size);
theta = Inf;

real_run = true;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backwards_error(A,b,y,theta)/norm(A,'fro')];
else
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);0];
end

%% LSQR

[~,lsqr]=sketch_and_precondition(A,b,20*n,trials,summary,true);

%% Iterative Sketching

[~,itsk]=iterative_sketching(A,b,20*n,trials,summary,true);

%% QR

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:trials, lsqr(:,j), 'LineWidth', 4, 'Color', "#0072BD"); hold on
    semilogy(0:trials, itsk(:,j), 'LineWidth', 4, 'Color', "#EDB120");
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        legend('S\&P ',...
            'IS',...
            'QR')
	if real_run
            saveas(gcf,'figs/sketch_precondition_forward.fig')
            saveas(gcf,'figs/sketch_precondition_forward.png')
	end
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $r$}(\mbox{\boldmath $x$})-\mbox{\boldmath $r$}(\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
        legend('S\&P ',...
            'IS',...
            'QR')
	if real_run
            saveas(gcf,'figs/sketch_precondition_residual.fig')
            saveas(gcf,'figs/sketch_precondition_residual.png')
	end
    elseif j == 3
        ylabel('Backward error $\mbox{BE}(\mbox{\boldmath $\widehat{x}$}_i)$')
        legend('S\&P ',...
            'IS',...
            'QR')
	if real_run
            saveas(gcf,'figs/sketch_precondition_backward.fig')
            saveas(gcf,'figs/sketch_precondition_backward.png')
	end
    end
end

%% Save

if real_run
    save('data/results_variants.mat', 'lsqr', 'itsk', 'qr_vals', 'trials')
end