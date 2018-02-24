%KEPLER Solves Kepler's equation for eccentric anomaly.

%   E = kepler(e,M) solves M = E - e*sin(E) for E

%   where E is eccentric anomaly in radians, e is eccentricity, 

%   and M is mean anomaly in radians.



function E = kepler(e,M)



i = 1;

tol = 1E-12;



B = cos(e) - ((pi/2) - e)*sin(e);



E(1) = M + (e*sin(M))/(B + M*sin(e));



while abs(E(i) - M - e*sin(E(i))) > tol   

   fE(i) = E(i) - e*sin(E(i)) - M;

   dfE(i) = 1 - e*cos(E(i));

   d2fE(i) = e*sin(E(i));

   

   A = 2*sqrt(abs((4*(dfE(i)^2)) - 5*fE(i)*d2fE(i)));

   

   E1(i+1) = E(i) - (5*fE(i)/(dfE(i) + A));

   E2(i+1) = E(i) - (5*fE(i)/(dfE(i) - A));

   

   if (abs(E1(i+1) - E(i))) < (abs(E2(i+1) - E(i)))

      E(i+1) = E1(i+1);

   else      

      E(i+1) = E2(i+1);

   end

   i = i + 1;

   if i > 100

      fprintf('\nKepler`s Eqn. is NOT converging!\n')

      return

   end   

end



E = E(i);