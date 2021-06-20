function Bke= BKsiEtaQ8(ksi,eta)
Bke (1,:)=[ (1-eta)*(2*ksi+eta)/4, -(1-eta)*ksi, (1-eta)*(2*ksi-eta)/4, (1-eta^2)/2, (1+eta)*(2*ksi+eta)/4, -(1+eta)*ksi, (1+eta)*(2*ksi-eta)/4, -(1-eta^2)/2];
Bke (2,:)=[ (1-ksi)*(ksi+2*eta)/4, -(1-ksi^2)/2, -(1+ksi)*(ksi-2*eta)/4, -(1+ksi)*eta, (1+ksi)*(ksi+2*eta)/4, (1-ksi^2)/2, -(1-ksi)*(ksi-2*eta)/4, -(1-ksi)*eta];
end