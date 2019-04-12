% user inputs params
re=100;
ny = 64; % number of elements in y
nx = ny ;  % number of elements in x
r = 1;  % aspect ratio
dx = 1/ny; % x spacing
dy = 1/nx;  % y spacing
% initializaation
[u, v, wMat,psiMat, uChng, vChng, wChng, psiChng] = deal(zeros(ny+1, nx+1));

% Boundary condition for wMat, top layer
for i=1:ny+1
    u(i,nx+1)=1;
    wMat(i,nx+1)=-2/dy;
end

uChng(1,nx+1)=1; % to enter the loop first time
while(max(max(abs(uChng)))>10^-6 || max(max(abs(vChng)))>10^-6 || max(max(abs(wChng)))>10^-6 || max(max(abs(psiChng)))>10^-6)
   wMatOld=wMat;
   psiMatOld=psiMat;
   uOld=u;
   vOld=v;

% second order upwind schemem for all nodes except the once on last two edges
for i=3:ny-1
    for j=3:nx-1
        if u(i,j)>=0 && v(i,j)>0
            k=u(i,j)/2/dx*(4*wMat(i-1,j)-wMat(i-2,j))+v(i,j)/2/dy*(4*wMat(i,j-1)-wMat(i,j-2));
            l=3/2*(u(i,j)/dx+v(i,j)/dy);

        elseif u(i,j)>=0 && v(i,j)<0
            k=u(i,j)/2/dx*(4*wMat(i-1,j)-wMat(i-2,j))+v(i,j)/2/dy*(wMat(i,j+2)-4*wMat(i,j+1));
            l=3/2*(u(i,j)/dx-v(i,j)/dy);

        elseif u(i,j)<0 && v(i,j)>=0
            k=u(i,j)/2/dx*(wMat(i+2,j)-4*wMat(i+1,j))+v(i,j)/2/dy*(4*wMat(i,j-1)-wMat(i,j-2));
            l=3/2*(-u(i,j)/dx+v(i,j)/dy);

        elseif u(i,j)<0 && v(i,j)<0
            k=u(i,j)/2/dx*(wMat(i+2,j)-4*wMat(i+1,j))+v(i,j)/2/dy*(wMat(i,j+2)-4*wMat(i,j+1));
            l=3/2*(-u(i,j)/dx-v(i,j)/dy);

        elseif u(i,j)>0 && v(i,j)==0
            k=u(i,j)/2/dx*(4*wMat(i-1,j)-wMat(i-2,j));
            l=3/2*u(i,j)/dx;
        end
        wMat(i,j)=(1/re*((wMat(i+1,j)+wMat(i-1,j))/dx(1,1)^2+1/r^2*(wMat(i,j+1)+wMat(i,j-1))/dy(1,1)^2)+k)/(l+2/re*(1/dx(1,1)^2+1/r^2/dy(1,1)^2));
    end
end

% First order upwind scheme for nodes lying on one layer internal layer
for i=2:ny
    for j=2:nx
        if i==2 || i==ny || j==2 || j==nx
             if u(i,j)>=0 && v(i,j)>0
               k=u(i,j)*wMat(i-1,j)/dx+v(i,j)*wMat(i,j-1)/dy;
               l=u(i,j)/dx+v(i,j)/dy;

             elseif u(i,j)>=0 && v(i,j)<0
               k=u(i,j)*wMat(i-1,j)/dx-v(i,j)*wMat(i,j+1)/dy;
               l=u(i,j)/dx-v(i,j)/dy;

             elseif u(i,j)<0 && v(i,j)>=0
               k=-u(i,j)*wMat(i+1,j)/dx+v(i,j)*wMat(i,j-1)/dy;
               l=-u(i,j)/dx+v(i,j)/dy;

             elseif u(i,j)<0 && v(i,j)<0
               k=-u(i,j)*wMat(i+1,j)/dx-v(i,j)*wMat(i,j+1)/dy;
               l=-u(i,j)/dx-v(i,j)/dy;

             elseif u(i,j)>0 && v(i,j)==0
               k=u(i,j)*wMat(i-1,j)/dx;
               l=u(i,j)/dx;
        end
        wMat(i,j)=(1/re*((wMat(i+1,j)+wMat(i-1,j))/dx^2+1/r^2*((wMat(i,j+1)+wMat(i,j-1))/dy^2))+k)/(l+2/re*(1/dx^2+1/r^2/dy^2));
        end
    end
end

% calculating stream function psi
for i=2:ny
    for j=2:nx
        psiMat(i,j)=(r^2*(psiMat(i+1,j)+psiMat(i-1,j))/dx^2+(psiMat(i,j+1)+psiMat(i,j-1))/dy^2+wMat(i,j))/(2*r^2/dx^2+2/dy^2);
    end
end

% calculating the u,v
for i=2:ny
    for j=2:nx
        u(i,j)=(psiMat(i,j+1)-psiMat(i,j-1))/(2*dy);
        v(i,j)=-(psiMat(i+1,j)-psiMat(i-1,j))/(2*dx);
    end
end

% Aaplying BCs to wMat on edges
% top layer
for i=1:ny+1
    wMat(i,nx+1)=-2*psiMat(i,nx)/dy^2-2/dy;
end
% right layer
for  j=2:nx
    wMat(ny+1,j)=r^2*(-2*psiMat(ny,j)/dx^2);
end
% bottom layer
for  i=2:ny
    wMat(i,1)=-2*psiMat(i,2)/dy^2;
end
% left layer
for j=2:nx
    wMat(1,j)=r^2*(-2*psiMat(2,j)/dx^2);
end

% calculating the changes in old and new matrics
uChng=max(max(abs(u-uOld)));
vChng=max(max(abs(v-vOld)));
wChng=max(max(abs(wMat-wMatOld)));
psiChng=max(max(abs(psiMat-psiMatOld)));
end

% linear vector with equal spacing
X1 = linspace(0,1,nx+1);
Y1 = linspace(0,1,ny+1);

vD2=zeros(nx+1,1);
for i=1:ny+1
    % v at x= D/2
    vD2(i,1)=v(i,nx/2+1);
end

uH2=zeros(nx+1,1);
for j=1:nx+1
    % u at y= H/2
    uH2(j,1)=u(ny/2+1,j);
end

% standard values
vStd = zeros(17,6);
vStd(:,1)=[0.0000;-.05906;-.07391;-.08864;-.10313;-.16914;-.22445;-.24533;.05454;.17527;.17507;.16077;.12317;.10890;.10091;.09233;.00000];
vStd(:,2)=[0.0000;-.12146;-.15663;-.19254;-.22847;-.23827;-.44993;-.38598;.05186;.30174;.30203;.28124;.22965;.20920;.19713;.18360;.00000];
vStd(:,3)=[0.0000;-.21388;-.27669;-.33714;-.39188;-.51550;-.42665;-.31966;.02526;.32235;.33075;.37095;.32627;.30353;.29012;.27485;.00000];
vStd(:,4)=[0.0000;-.39017;-.47425;-.52357;-.54053;-.44307;-.37401;-.31184;.00999;.28188;.29030;.37119;.42768;.41906;.40917;.39560;.00000];
vStd(:,5)=[0.00000;-.49774;-.55069;-.55408;-.52876;-.41442;-.36214;-.30018;.00945;.27280;.28066;.35368;.42951;.43648;.43329;.42447;.0000];
vStd(:,6)=[0.0000;-.53858;-.55216;-.52347;-.48590;-.41050;-.36213;-.30448;.00824;.27348;.28117;.35060;.41824;.43564;.44030;.43979;.0000];
xStd=[1.0000;.9688;.9609;.9531;.9453;.9063;.8594;.8047;.5000;.2344;.2266;.1563;.0938;.0781;.0703;.0652;.0000];

uStd = zeros(17,6);
uStd(:,1)=[1.0000;.84123;.78871;.73722;.68717;.23151;.00332;-.13641;-.20581;-.21090;-.15662;-.10150;-.06434;-.04775;-.04192;-.03717;.00000];
uStd(:,2)=[1.0000;.75837;.68439;.61756;.55892;.29093;.16256;.02135;-.11477;-.17119;-.32726;-.24299;-.14612;-.10338;-.09266;-.08186;.00000];
uStd(:,3)=[1;0.65928;0.57492;.51117;.46604;.33304;.18719;.05702;-.0608;-.10648;-.27805;-.38289;-.2973;-.2233;-.20196;-.18109;0];
uStd(:,4)=[1;0.53236;0.48296;0.46547;0.46101;0.34682;0.19791;0.07156;-0.04272;-0.86636;-0.24427;-0.34323;-0.41933;-0.37827;-0.35344;-0.32407;0];
uStd(:,5)=[1;0.48223;0.4612;0.45992;0.46036;0.33556;0.20087;0.08183;-0.03039;-0.07404;-0.22855;-0.33050;-0.40435;-0.43643;-0.42901;-0.41165;0];
uStd(:,6)=[1;0.47244;0.47048;0.47323;0.47167;0.34228;0.20591;0.08342;-0.038;-0.07503;-0.23176;-0.32393;-0.38324;-0.43025;-0.43590;-0.43154;0];
yStd=[1.0000;.9766;.9688;.9609;.9531;.8516;.7344;.6172;.5000;.4531;.2813;.1719;.1016;.0703;.0625;.0547;.0000];

plot(Y1,uH2)
hold on

plot(yStd,uStd(:,1),'*')
ylabel('u/u')
xlabel('y/H')

plot(X1,vD2)
hold on

plot(xStd,vStd(:,1),'*')
ylabel('v/v')
xlabel('x/D')    
