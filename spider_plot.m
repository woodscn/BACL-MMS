function out = spider_plot()
dat = load('random_heat_exact_best.dat');
dat(:,2:3) = dat(:,2:3)*pi/180.;
dat(:,4:5) = log10(abs(dat(:,4:5)));
dat(:,8:9) = log10(abs(dat(:,8:9)));
labels = {'Problem #','\theta (rad)','\phi (rad)',...
    'log(A)','log(B)','C','\mu','log(\alpha)',...
    'log(Residual)'};
figure(1)
parallelcoords(dat,'Labels',labels)
ylim([-10,10])
datE = load('Eulerdats.mat');
datE = datE.datE;
for n = 1:57
    norm_facs = normalizing_factors(datE(n,1)+1);
    datE(n,3:end) = log10(...
        abs(datE(n,3:end))./norm_facs...
        +(1e-19));
end
labelsE = {'Problem #','\theta (rad)','log(R(rho))',...
    'log(R(rho u))', 'log(R(rho v))', 'log(R(rho e))'};
figure(2)
parallelcoords(datE,'Labels',labelsE)
ylim([-10 5])
out = 0;
end
function out = primtocons(in)
p = in(1); d = in(2); u = in(3); v = in(4); w = in(5);
out = [d,d*u,d*v,d*w,d*(.5*(u^+v^2+w^2)+p/(.4*d))];
end
function out = states_func()
out(1,:,:) = [[1,1,0,0,0];[.30313,.42632,.92745,0,0];...
    [.30313,.26557,.92745,0,0];[.1,.125,0,0,0]];
out(2,:,:) = [[.4,1,-2,0,0];[.00189,.02185,0,0,0];...
    [.00189,.02185,0,0,0];[.4,1,2,0,0]];
out(3,:,:) = [[1000,1,0,0,0];[460.894,.57506,19.5975,0,0];...
    [460.894,5.99242,19.5975,0,0];[.01,1,0,0,0]];
out(4,:,:) = [[.01,1,0,0,0];[46.0950,5.99242,-6.19633,0,0];...
    [46.0950,.57511,-6.19633,0,0];[100,1,0,0,0]];
end
function out = normalizing_factors(n)
max_cons = [0,0,0,0,0];
for ind = 1:5
    states = states_func();
    for inda = 1:4
        cons_states(inda,:) = primtocons(reshape(states(n,inda,:),[1,5]));
    end
    max_cons(ind) = max(abs(cons_states(:,ind)));
end
for ind = 1:5 
    if max_cons(ind) == 0
        max_cons(ind) = 1;
    end
end
out = max_cons;
out = [out(1),out(2),out(3),out(5)];
end
