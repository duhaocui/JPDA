function xk1=processmodel_2body(dt,MU,tk,xk)
@(t,x)twoBody(t,x);


opt = odeset('reltol',1e-6,'abstol',1e-6);
% [~,xx]=ode45(@twoBody,[tk,tk+dt],xk,opt);
% xk1=xx(end,:)';
[ r, v, Ehat ] = FnG(tk, tk+dt, xk(1:3), xk(4:6), MU);


xk1=[r(:);v(:)];
end