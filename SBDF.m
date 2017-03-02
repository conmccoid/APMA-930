%% Steady Navier Stokes
%
% Solve the driven cavity problem for Stokes flow in a wedge
% using the streamfunction/vorticity formulation.
%
% Builds a big system matrix
%
    close all;  clear; 
    clc;
    
%
% Problem parameters:
    U = -1;
    Rmax = 1;
    alpha = pi/2;
    Re = 1000;
    dt = 1e-2;
    a = 2*dt/(3*Re);

%
% Set up finite difference grid
    M =100; 
    dr = Rmax/(M-1);
    N =100; 
    dth = alpha/(N-1);
    [rg, thg] = meshgrid(0: dr :Rmax, ...
                         alpha: -dth: 0);                     
%
% Unknowns and numbering
    numUn = M*N;
    nP = reshape(1:numUn, size(rg));
    nO = reshape(numUn+1:2*numUn, size(rg));
    numUn = 2*numUn;
%
% Build system matrices and rhs
%     PsiOmSys = SystemMat(numUn, nP, nO, M, N, alpha, dr, dth);
    PsiOmSys = SystemMatCompSBDF(numUn, nP, nO, M, N, alpha, dr, dth, a);
    PsiOmSys = sparse(PsiOmSys);
    A = PsiOmSys;
    [LL, UU, PP, QQ, RR] = lu(A);

    spparms('spumoni', 0)  % write out information about sparse algorithms
%     spy(PsiOmSys)
%     drawnow

%%    

% Initialization
    tol = 1e-3;

    psivort_old = zeros(numUn,1);
    psivort = zeros(numUn,1);
    J_old = zeros(numUn,1);
    
    TestPsi = 1; TestOm = 1;
    IterMax = 1000; Iter = 0;

    t = 0;
    T = 2000*dt;
%
% Solve
% while Iter<IterMax && (TestPsi>tol || TestOm>tol)
while t<T && isnan(TestOm)~=1
    t = t+dt;
%     tic
    rhs_new = ConstructRhsSBDF(numUn, nP, nO, M, N, Rmax, dr, dth, dt, U, psivort, psivort_old);
    c = PP * (RR \ rhs_new);
    psivort_new = QQ * ( UU \ ( LL \ c ) );
%     t = toc;
%     disp(['Time taken for linear system solve = ', num2str(t)]);

    TestPsi = norm(psivort_new(1:numUn/2) - psivort(1:numUn/2))/norm(psivort_new(1:numUn/2));
    TestOm = norm(psivort_new((numUn/2+1):end) - psivort((numUn/2+1):end))/norm(psivort_new((numUn/2+1):end));
    Iter = Iter + 1;
    disp(['Iteration: ',num2str(Iter),' TestPsi: ',num2str(TestPsi),' TestOm: ',num2str(TestOm)]);
    
    psivort_old = psivort;
    psivort = psivort_new;
    
end

    psi = reshape(psivort(1:numUn/2), size(rg));
    omega = reshape(psivort(numUn/2+1:numUn), size(rg)); 
  
%%
% Plot

    figure()
    subplot(1, 2, 1)
        pcolor(rg.*cos(thg), rg.*sin(thg), psi); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([-Rmax Rmax -Rmax Rmax])
        axis square
    subplot(1, 2, 2)
        pcolor(rg.*cos(thg), rg.*sin(thg), omega); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity')
        axis([0 Rmax 0 Rmax])
        axis square
%%
%  Look for Eddies!
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 40); 
%         c = linspace(0, 4.4d-4, 40);
%         contour(rg.*cos(thg), rg.*sin(thg), psi, c); 
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 Rmax 0 Rmax])
        axis square
    subplot(1, 2, 2)
        c = linspace(-20, 20, 40);
        contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        shading flat;  colormap(jet);  
        xlabel('r')
        ylabel('\theta')
        title('Vorticity')
        axis([0 Rmax 0 Rmax])
        axis square