function Bke= BKsiEtaQ9(ksi,eta)
Bke (1,:)=[ (1-eta)*(1-2*ksi)*eta/4, (1-eta)*ksi*eta, -(1-eta)*(1+2*ksi)*eta/4, (1+2*ksi)*(1-eta^2)/2, (1+eta)*(1+2*ksi)*eta/4, -(1+eta)*ksi*eta, -(1+eta)*(1-2*ksi)*eta/4, -(1-eta^2)*(1-2*ksi)/2, -2*(1-eta^2)*ksi];
Bke (2,:)=[ (1-ksi)*(1-2*eta)*ksi/4, -(1-ksi^2)*(1-2*eta)/2, -(1+ksi)*(1-2*eta)*ksi/4, -(1+ksi)*ksi*eta, (1+ksi)*(1+2*eta)*ksi/4, (1-ksi^2)*(1+2*eta)/2, -(1-ksi)*(1+2*eta)*ksi/4, (1-ksi)*eta*ksi, -2*(1-ksi^2)*eta];
end