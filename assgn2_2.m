%Runs through StokesEddies a number of times to test effect of M and N on
%accuracy of physical boundary condition

j = 1; 
% errM = zeros(length(50:10:200),1);
err2 = zeros(length(50:10:200),1); 
err1 = err2; 
errinf = err2;
% N = 50;
for N = 50:10:200
    M = N;
    StokesEddies
    diff = zeros(M,1);
    diff([M-2;M-1;M]) = [1/(2*dr);-4/(2*dr);3/(2*dr)];
%     errM(j) = norm(psi*diff - 1)/norm(ones(M,1));
    err2(j) = norm(psi*diff - 1)/norm(ones(M,1));    
    err1(j) = norm(psi*diff - 1,1)/norm(ones(M,1),1);    
    errinf(j) = norm(psi*diff - 1,inf)/norm(ones(M,1),inf);
    j = j+1;
end

loglog(50:10:200,err2,'ko--',50:10:200,err1,'ks--',50:10:200,errinf,'k^--','MarkerSize',10,'Linewidth',5)
xlabel('M and N')
ylabel('Error')
legend('2-norm','1-norm','Inf-norm')
set(gca,'FontSize',26,'Linewidth',5)
axis([50,200,1e-2,1])