function Xsat0=getorbits(n)


% [a0,e0,i0,omg0,Omg0,M0]=cart2orb(7000,0,0,0,-1.0374090357,7.4771288355);
mue = 398601.2;

a0=[7000,8000];
e0=[0.1,0.4];
M0=[0,15]*pi/180;
w0=[10,85]*pi/180;
i0=[-40,40]*pi/180;
Omega0=[0,80]*pi/180;

i=1;
Xsat0=zeros(n,6);
while(i<=n)
    a = a0(1) + (a0(2)-a0(1)).*rand(1,1);
    ecc =  e0(1) + (e0(2)-e0(1)).*rand(1,1);
    if ecc<0
        ecc
        continue
    end
    if a*(1-ecc)<6700
        continue
    end
    inc= i0(1) + (i0(2)-i0(1)).*rand(1,1);
    Omega= Omega0(1) + (Omega0(2)-Omega0(1)).*rand(1,1);
    w=w0(1) + (w0(2)-w0(1)).*rand(1,1);
    M=M0(1) + (M0(2)-M0(1)).*rand(1,1);
    if M<0
%         M
        continue
    end
    
    E=kepler(ecc,M);
%     E
    XX= OE2XYZ([a, ecc, E, w, inc, Omega,mue]);
    r=XX(1:3);
    v=XX(4:6);
    
%     if a*(1-ecc)>6500
        Xsat0(i,:)=[r',v'];
        i=i+1;
%     end

end

end