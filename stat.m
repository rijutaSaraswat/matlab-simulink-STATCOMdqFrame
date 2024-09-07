clear all
clc

f0=60; Wb=2*pi*f0; W0=Wb;
Rs=0.01; Xs=0.15; rp=78.7;
gp=1/rp; bc=1.136;
k=2*sqrt(6)/pi;

%STATCOM initial values
vdc=0.7; iQ=0; iD=1; vD=0; vQ=1;

% disp('Enter iRref=1 for inductive else iRref=-1 for capacitive');
% iRref=input('Enter iRref value?\n');
% disp('Enter 0 for PI controller else 1 for PI controller with non-linear feedback');
% Controller=input('Enter Controller value?\n');
iRref=1;
Controller=1;
teta=atan(vD/vQ);

%TYPE-2 SUB-OPTIMAL CONTROLLER PARAMETERS
g=2;
Tw=0.01;
kp=0.33;
ki=10*kp;

%TYPE-2 CONTROLLER PARAMETERS
g=2.54;
Tw=0.018;
kp=0.69;
ki=16.45;

althetadeg=-0.5;
altheta=(pi/180)*althetadeg;

%INITIAL VALUE PROGRAM
% unknowns iD; iQ;  vdc;  altheta;
func(1,1)=Wb/Xs*(-Rs*iD-((W0/Wb)*Xs*iQ)+vD-(k*vdc*sin(altheta)));
func(2,1)=Wb/Xs*(-Rs*iQ+((W0/Wb)*Xs*iD)+vQ-(k*vdc*cos(altheta)));
func(3,1)=Wb/bc*(-gp*vdc+(k*iD*sin(altheta))+(k*iQ*cos(altheta)));
func(4,1)=-iD*cos(teta)+iQ*sin(teta)-iRref;

while  max(abs(func))>1e-12 
func(1,1)=Wb/Xs*(-Rs*iD-((W0/Wb)*Xs*iQ)+vD-(k*vdc*sin(altheta)));
func(2,1)=Wb/Xs*(-Rs*iQ+((W0/Wb)*Xs*iD)+vQ-(k*vdc*cos(altheta)));
func(3,1)=Wb/bc*(-gp*vdc+(k*iD*sin(altheta))+(k*iQ*cos(altheta)));
func(4,1)=-iD*cos(teta)+iQ*sin(teta)-iRref;

% Jacobian matrix elements
jacob(1,1)=-(Wb/Xs*Rs);
jacob(1,2)=-W0;
jacob(1,3)=-(Wb/Xs)*k*sin(altheta);
jacob(1,4)=(Wb/Xs)*(-k*vdc*cos(altheta));

jacob(2,1)=W0;
jacob(2,2)=-(Wb/Xs*Rs);
jacob(2,3)=-(Wb/Xs)*k*cos(altheta);
jacob(2,4)=(Wb/Xs)*(k*vdc*sin(altheta));

jacob(3,1)=(Wb/bc)*k*sin(altheta);
jacob(3,2)=(Wb/bc)*k*cos(altheta);
jacob(3,3)=-gp*(Wb/bc);
jacob(3,4)=(Wb/bc)*(k*iD*cos(altheta)-k*iQ*sin(altheta));

jacob(4,1)=-cos(teta);
jacob(4,2)=sin(teta);
jacob(4,3)=0;
jacob(4,4)=0;
 
deltaX=-jacob^(-1)*func;
    iD=iD+deltaX(1,1);
    iQ=iQ+deltaX(2,1);
    vdc=vdc+deltaX(3,1);
    altheta=altheta+deltaX(4,1);
end

% initial conditions for the TYPE-2 DQ-STATCOM
iD;
iQ;
vdc;
altheta;
alfa=altheta-teta;
alfa
alfadeg=alfa*(180/pi);
alfadeg
iR=-iD*cos(teta)+iQ*sin(teta)
iP=iD*sin(teta)+iQ*cos(teta)

[sizes x0 xstored]=statcom_model_Type2([],[],[],0);
[a,b,c,d]=linmodv5('statcom_model_Type2',x0);
eigval=eig(a)
