function fid = fidelityRho(rho1, rho2)

rho1 = rho1/trace(rho1);
rho2 = rho2/trace(rho2);



fid = trace(sqrtm(sqrtm(rho2)*rho1*sqrtm(rho2)));
