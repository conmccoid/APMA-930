function [r] = ConstructRhsSBDF(numUn, nP, nO, M, N, R, dr, dth, dt, U, psivort, psivort_old)
    r = zeros(numUn, 1);
%
%  Boundary conditions at right: omega given
    for jrow = 2: N-1
        ijO = nO(jrow, M);
        r(ijO) = U/R + 3*U/dr;
    end  
    
%  Vorticity
    for jrow = 2:N-1
        for icol = 2:M-1
            r_i = (icol-1)*dr;
            ijO = nO(jrow,icol);
            ijpO = nO(jrow+1,icol);
            ijmO = nO(jrow-1,icol);
            ipjO = nO(jrow,icol+1);
            imjO = nO(jrow,icol-1);
            ijP = nP(jrow,icol);
            ijpP = nP(jrow+1,icol);
            ijmP = nP(jrow-1,icol);
            ipjP = nP(jrow,icol+1);
            imjP = nP(jrow,icol-1);
            J_new = ((psivort(ijpO)-psivort(ijmO))*(psivort(ipjP)-psivort(imjP)) - ...
                (psivort(ipjO)-psivort(imjO))*(psivort(ijpP)-psivort(ijmP)))/(4*dr*dth*r_i);
            J_old = ((psivort_old(ijpO)-psivort_old(ijmO))*(psivort_old(ipjP)-psivort_old(imjP)) - ...
                (psivort_old(ipjO)-psivort_old(imjO))*(psivort_old(ijpP)-psivort_old(ijmP)))/(4*dr*dth*r_i);
            r(ijO) = (2*dt/3)*(J_old - 2*J_new);
            r(ijO) = r(ijO)+(4*psivort(ijO) - psivort_old(ijO))/3;
        end
    end
end