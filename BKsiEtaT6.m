function Bke= BKsiEtaT6(ksi,eta)
lambda=1-ksi-eta;
Bke (1,:)=[ 1-4*lambda, 4*(lambda-ksi), -1+4*ksi, eta*4, 0, -4*eta];
Bke (2,:)=[ 1-4*lambda, -4*ksi, 0, ksi*4, -1+4*eta, 4*(lambda-eta)];
end
