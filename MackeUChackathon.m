op=NaN;
while op~=0 && op~=1 && op~=2 && op~=3 && op~=4
    op=input('Select an option:\n [0]:Exit\n [1]:Isen flow cal\n [2]:Wind Tunnel\n [3]:Cavity flow\n [4]:POD\nSelection:');
    if op==1
        op=NaN;
%         disp('Option1')
clear variables
        g=input('Input Ratio of Specific Heats: ');
        M=input('Input Mach Number: ');
        m=(0:0.01:5);
        ToT=1+0.5*(g-1)*m.^2;
        tot=1+0.5*(g-1)*M^2;
        PoP=ToT.^(g/(g-1));
        pop=tot^(g/(g-1));
        RoR=ToT.^(1/(g-1));
        ror=tot^(1/(g-1));
        figure(1)
        subplot(1,3,1)
        plot(m,ToT)
        hold on
        scatter(M,tot)
        hold off
        title(['M=',num2str(M),', ToT=',num2str(tot)])
        xlabel('Mach Number')
        ylabel('To/T')
        subplot(1,3,2)
        plot(m,PoP)
        hold on
        scatter(M,pop)
        hold off
        title(['M=',num2str(M),', PoP=',num2str(pop)])
        xlabel('Mach Number')
        ylabel('Po/P')
        subplot(1,3,3)
        plot(m,RoR)
        hold on
        scatter(M,ror)
        hold off
        title(['M=',num2str(M),', \rhoo/\rho=',num2str(ror)])
        xlabel('Mach Number')
        ylabel('\rhoo/\rho')
    elseif op==2
        op=NaN;
%         disp('Option2')
        clear variables
 
menu=4;
while menu~=1 && menu~=2 && menu~=3 
    menu=input('Select Simulation: \n  [1]: Simple Duct \n  [2]: Obstructed Duct \n  [3]: Turbine Cascade \nSelection:');
end
 
% DRIVING PARAMETERS %
U=1;%inlet velocity
Re=input('Reynolds number: ');%Reynolds number
gd=80;%grid density, [points/unit length]
dt=0.0001;% dt (time step)
inner=1;%iterations for vorticity 
sc=9e-6;%stop condition, guaranties steady state condition
q=1000;%print results every 'q' time steps (simulation resolution = dt*q)
dnorth=0.25;%diameters upper half
dsouth=0.25;%diameters lower half
dwest=0.5;%diameters upstream
deast=1;%diameters downstream
ly=dnorth+dsouth;%duct y-dimension
lx=dwest+1+deast;%duct x-dimension
px=lx*gd;%grid points in x direction
py=ly*gd+1;%grid points in y direction
X=zeros(py,px);%X grid matrix
Y=zeros(py,px);%Y grid matrix
 
% GRID GENERATION %
if menu==1%dense grid near boundaries
    x=linspace(0,lx,px);
    yc=ly/2;
    Y((py+1)/2,:)=yc;
    for j=1:py
        for i=1:px
           X(j,i)=x(1,i); 
        end
    end
    x=linspace(0,1,(py+1)/2);
    y=yc*x.^1.5;
    for i=1:px
       for j=1:(py-1)/2
          Y(j,i)=y(1,j);
          dy=yc-y(1,j);
          Y(py+1-j,i)=yc+dy;
       end
    end
elseif menu==2
    naca=2316;%NACA for case [2]
    Ndigits = dec2base(naca,10) - '0';
    if length(Ndigits)==3
        Ndigits=[0 Ndigits];
    elseif length (Ndigits)==2
        Ndigits=[0 0 Ndigits];
    elseif length(Ndigits)==1
        Ndigits=[0 0 0 Ndigits];
    end
    m=Ndigits(1,1)/100;%maximum camber
    p=Ndigits(1,2)/10;%location of maximum camber
    t=((Ndigits(1,3)*10)+Ndigits(1,4))/100;%thickness 
    beta=linspace(0,pi,gd);
    yt=zeros(1,gd);
    yc=yt;xu=yt;xl=yt;yu=yt;yl=yt;theta=yt;
    xc=linspace(0,1,gd);
    %airfoil
    for i=1:gd
        yt(1,i)=(t/0.2)*(0.2969*(xc(1,i)^0.5)-0.1260*xc(1,i)-0.3516*(xc(1,i)^2)+0.2843*(xc(1,i)^3)-0.1015*(xc(1,i)^4));
        if xc(1,i)<p
            yc(1,i)=(m/p^2)*(2*p*xc(1,i)-xc(1,i)^2);
            theta(1,i)=atan((2*m/(p^2))*(p-xc(1,i)));
        elseif xc(1,i)>=p
            yc(1,i)=(m/(1-p)^2)*(1-2*p+2*p*xc(1,i)-xc(1,i)^2);
            theta(1,i)=atan((2*m/(1-p)^2)*(p-xc(1,i)));
        end
        xu(1,i)=dwest+xc(1,i)-yt(1,i)*sin(theta(1,i));
        xl(1,i)=dwest+xc(1,i)+yt(1,i)*sin(theta(1,i));
        yu(1,i)=yc(1,i)+yt(1,i)*cos(theta(1,i));
        yl(1,i)=ly+yc(1,i)-yt(1,i)*cos(theta(1,i));
    end
    yc=ly/2;
    Y((py+1)/2,:)=yc;
    xwest=linspace(0,dwest,(dwest*gd)+1);
    xeast=linspace(dwest+1,lx,(deast*gd)+1);
    xnorth=[xwest(1,1:dwest*gd),xl,xeast(1,2:(deast*gd)+1)];
    xsouth=[xwest(1,1:dwest*gd),xu,xeast(1,2:(deast*gd)+1)];
    xc=linspace(0,lx,lx*gd);
    X((py+1)/2,:)=xc;
    for j=1:py
       if j<(py+1)/2
          for i=1:px
              X(j,i)=xsouth(1,i);
          end
       elseif j>(py+1)/2
           for i=1:px
               X(j,i)=xnorth(1,i);
           end
       end
    end
    for i=1:px
       if i<(dwest*gd)+1
          ynorth=linspace(ly/2,ly,(dnorth*gd)+1);
          ysouth=linspace(0,ly/2,(dsouth*gd)+1);
          for j=1:(py-1)/2
              Y(((py+1)/2)+j,i)=ynorth(1,j+1);
              Y(j,i)=ysouth(1,j);
          end
       elseif i>dwest*gd && i<((dwest+1)*gd)+1
          ynorth=linspace(ly/2,ly,(dnorth*gd)+1);
          ysouth=linspace(yu(1,i-(dwest*gd)),ly/2,(dsouth*gd)+1);
          for j=1:(py-1)/2
              Y(((py+1)/2)+j,i)=ynorth(1,j+1);
              Y(j,i)=ysouth(1,j);%Y(j,2*dwest*gd+gd+1-i)=ysouth(1,j);%
          end
       elseif i>=((dwest+1)*gd)+1
          ynorth=linspace(ly/2,ly,(dnorth*gd)+1);
          ysouth=linspace(yu(1,gd),ly/2,(dsouth*gd)+1);
          for j=1:(py-1)/2
              Y(((py+1)/2)+j,i)=ynorth(1,j+1);
              Y(j,i)=ysouth(1,j);
          end
       end
    end
    for i=1:size(Y,2)
       ys=linspace(Y(1,i),ly,size(Y,1));
       for j=1:size(Y,1)
          Y(j,i)=ys(1,j); 
       end
    end
elseif menu==3
    naca=3312;%NACA for case [3]
    dnorth=0.25;%diameters upper half
    dsouth=0.25;%diameters lower half
    dwest=0;%diameters upstream
    deast=1;%diameters downstream
    ly=dnorth+dsouth;%duct y-dimension
    lx=dwest+1+deast;%duct x-dimension
    px=lx*gd;%grid points in x direction
    py=ly*gd+1;%grid points in y direction
    X=zeros(py,px);%X grid matrix
    Y=zeros(py,px);%Y grid matrix
    Ndigits = dec2base(naca,10) - '0';
    if length(Ndigits)==3
        Ndigits=[0 Ndigits];
    elseif length (Ndigits)==2
        Ndigits=[0 0 Ndigits];
    elseif length(Ndigits)==1
        Ndigits=[0 0 0 Ndigits];
    end
    m=Ndigits(1,1)/100;%maximum camber
    p=Ndigits(1,2)/10;%location of maximum camber
    t=((Ndigits(1,3)*10)+Ndigits(1,4))/100;%thickness 
    beta=linspace(0,pi,gd);
    yt=zeros(1,gd);
    yc=yt;xu=yt;xl=yt;yu=yt;yl=yt;theta=yt;
    xc=linspace(0,1,gd);
    %airfoil
    for i=1:gd
        yt(1,i)=(t/0.2)*(0.2969*(xc(1,i)^0.5)-0.1260*xc(1,i)-0.3516*(xc(1,i)^2)+0.2843*(xc(1,i)^3)-0.1015*(xc(1,i)^4));
        if xc(1,i)<p
            yc(1,i)=(m/p^2)*(2*p*xc(1,i)-xc(1,i)^2);
            theta(1,i)=atan((2*m/(p^2))*(p-xc(1,i)));
        elseif xc(1,i)>=p
            yc(1,i)=(m/(1-p)^2)*(1-2*p+2*p*xc(1,i)-xc(1,i)^2);
            theta(1,i)=atan((2*m/(1-p)^2)*(p-xc(1,i)));
        end
        xu(1,i)=dwest+xc(1,i)-yt(1,i)*sin(theta(1,i));
        xl(1,i)=dwest+xc(1,i)+yt(1,i)*sin(theta(1,i));
        yu(1,i)=yc(1,i)+yt(1,i)*cos(theta(1,i));
        yl(1,i)=ly+yc(1,i)-yt(1,i)*cos(theta(1,i));
    end
    xl=xu;%%*
    yc=ly/2;
    Y((py+1)/2,:)=yc;
    xwest=linspace(0,dwest,(dwest*gd)+1);
    xeast=linspace(dwest+1,lx,(deast*gd)+1);
    xnorth=[xwest(1,1:dwest*gd),xl,xeast(1,2:(deast*gd)+1)];
    xsouth=[xwest(1,1:dwest*gd),xu,xeast(1,2:(deast*gd)+1)];
    xc=linspace(0,lx,lx*gd);
    X((py+1)/2,:)=xc;
    for j=1:py
       if j<(py+1)/2
          for i=1:px
              X(j,i)=xsouth(1,i);
          end
       elseif j>(py+1)/2
           for i=1:px
               X(j,i)=xnorth(1,i);
           end
       end
    end
    for i=1:px
       if i<(dwest*gd)+1
          ynorth=linspace(ly/2,ly,(dnorth*gd)+1);
          ysouth=linspace(0,ly/2,(dsouth*gd)+1);
          for j=1:(py-1)/2
              Y(((py+1)/2)+j,i)=ynorth(1,j+1);
              Y(j,i)=ysouth(1,j);
          end
       elseif i>dwest*gd && i<((dwest+1)*gd)+1
          ynorth=linspace(ly/2,yl(1,i-(dwest*gd)),(dnorth*gd)+1);
          ysouth=linspace(yu(1,i-(dwest*gd)),ly/2,(dsouth*gd)+1);
          for j=1:(py-1)/2
              Y(((py+1)/2)+j,i)=ynorth(1,j+1);
              Y(j,i)=ysouth(1,j);
          end
       elseif i>=((dwest+1)*gd)+1
          ynorth=linspace(ly/2,yl(1,gd),(dnorth*gd)+1);
          ysouth=linspace(yu(1,gd),ly/2,(dsouth*gd)+1);
          for j=1:(py-1)/2
              Y(((py+1)/2)+j,i)=ynorth(1,j+1);
              Y(j,i)=ysouth(1,j);
          end
       end
    end
end
pg=3;
while pg~=1 && pg~=0
    pg=input('\Plot grid? \n  [1]:Yes \n  [0]:No \n Selection: ');
    if pg==1
        %Plot Grid
        figure (1)
        for j=1:py
        plot(X(j,:),Y(j,:),'r')
        hold on
        end
        for i=1:px
        plot(X(:,i),Y(:,i),'r')
        end
        hold off
        title('Grid')
        xlabel('X')
        ylabel('Y')
        axis equal
        xlim([0 lx])
        ylim([0 ly])
    end
end
 
% METRIC COEFFICIENTS %
nz=px;%grid points in zeta direction
ne=py;%grid points in eta direction
J=zeros(ne,nz);%Jacobian field
I=zeros(ne,nz);%inverse Jacobian field
Zx=zeros(ne,nz);Zy=Zx;%Zeta_x and Zeta_y fields
Ex=zeros(ne,nz);Ey=Ex;%Eta_x and Eta_y fields
zeta=linspace(0,lx,nz);%zeta coordinate points
eta=linspace(0,ly,ne);%eta coordinate points
dz=zeta(1,2)-zeta(1,1);%dz (spacing in zeta direction)
de=eta(1,2)-eta(1,1);%de (spacing in eta direction)
%MC Computation
for j=1:ne
   for i=1:nz
       if i==1 && j==1
           xz=(X(j,i+1)-X(j,i))/(dz);
           ye=(Y(j+1,i)-Y(j,i))/(de);
           yz=(Y(j,i+1)-Y(j,i))/(dz);
           xe=(X(j+1,i)-X(j,i))/(de);
       elseif i==1 && j==ne
           xz=(X(j,i+1)-X(j,i))/(dz);
           ye=(Y(j,i)-Y(j-1,i))/(de);
           yz=(Y(j,i+1)-Y(j,i))/(dz);
           xe=(X(j,i)-X(j-1,i))/(de);
       elseif i==nz && j==1
           xz=(X(j,i)-X(j,i-1))/(dz);
           ye=(Y(j+1,i)-Y(j,i))/(de);
           yz=(Y(j,i)-Y(j,i-1))/(dz);
           xe=(X(j+1,i)-X(j,i))/(de);
       elseif i==nz && j==ne
           xz=(X(j,i)-X(j,i-1))/(dz);
           ye=(Y(j,i)-Y(j-1,i))/(de);
           yz=(Y(j,i)-Y(j,i-1))/(dz);
           xe=(X(j,i)-X(j-1,i))/(de);
       elseif i==1 && 1<j && j<ne
           xz=(X(j,i+1)-X(j,i))/(dz);
           ye=(Y(j+1,i)-Y(j-1,i))/(2*de);
           yz=(Y(j,i+1)-Y(j,i))/(dz);
           xe=(X(j+1,i)-X(j-1,i))/(2*de);
       elseif j==1 && 1<i && i<nz
           xz=(X(j,i+1)-X(j,i-1))/(2*dz);
           ye=(Y(j+1,i)-Y(j,i))/(de);
           yz=(Y(j,i+1)-Y(j,i-1))/(2*dz);
           xe=(X(j+1,i)-X(j,i))/(de);
       elseif i==nz && 1<j && j<ne
           xz=(X(j,i)-X(j,i-1))/(dz);
           ye=(Y(j+1,i)-Y(j-1,i))/(2*de);
           yz=(Y(j,i)-Y(j,i-1))/(dz);
           xe=(X(j+1,i)-X(j-1,i))/(2*de);
       elseif j==ne && 1<i && i<nz
           xz=(X(j,i+1)-X(j,i-1))/(2*dz);
           ye=(Y(j,i)-Y(j-1,i))/(de);
           yz=(Y(j,i+1)-Y(j,i-1))/(2*dz);
           xe=(X(j,i)-X(j-1,i))/(de);
       else
           xz=(X(j,i+1)-X(j,i-1))/(2*dz);
           ye=(Y(j+1,i)-Y(j-1,i))/(2*de);
           yz=(Y(j,i+1)-Y(j,i-1))/(2*dz);
           xe=(X(j+1,i)-X(j-1,i))/(2*de);
       end
       I(j,i)=xz*ye-yz*xe;
       Zx(j,i)=ye/I(j,i);
       Zy(j,i)=-xe/I(j,i);
       Ex(j,i)=-yz/I(j,i);
       Ey(j,i)=xz/I(j,i);
       J(j,i)=1/I(j,i);
   end
end
pg=3;
while pg~=1 && pg~=0
    pg=input('\Plot Jacobian? \n  [1]:Yes \n  [0]:No \n Selection: ');
    if pg==1
        %plot Jacobian
        figure(2)
        surf(X,Y,J,'EdgeColor','none');
        view(0,90);
        colormap jet;
        title('Jacobian');
        xlabel('X');
        ylabel('Y');
        colorbar
        axis equal
        xlim([0 lx])
        ylim([0 ly])
    end
end
 
% SOLVER %
T=0;
u=zeros(ne,nz);%u-velocity matrix
v=zeros(ne,nz);%v-velocity matrix
vort=u;%vorticity matrix
psi=zeros(2,ne,nz);%stream matrix
ppsi=zeros(ne,nz);%stream matrix (plot purposes)
omega=u;%Omega matrix
Vr=u;%resultant velocity matrix
%Boundary Conditions
%Left
u(:,1)=U;
v(:,1)=0;
for j=2:ne-1
   psi(1,j,1)=U*Y(j,1);
   psi(2,j,1)=psi(1,j,1);
   omega(j,1)=(1/Re)*(Zx(j,1)*((v(j,2)-v(j,1))/dz)-Zy(j,1)*((u(j,2)-u(j,1))/dz));
   vort(j,1)=Re*omega(j,1);
end
%Right
u(:,nz)=u(:,nz-1);
v(:,nz)=v(:,nz-1);
for j=2:ne-1
   psi(1,j,nz)=psi(1,j,nz-1);
   psi(2,j,nz)=psi(1,j,nz);
   omega(j,nz)=(1/Re)*(Ex(j,nz)*((v(j+1,nz)-v(j-1,nz))/(2*de))-Ey(j,nz)*((u(j+1,nz)-u(j-1,nz))/(2*de)));
   vort(j,nz)=Re*omega(j,nz);
end
%Top
u(ne,:)=0;
v(ne,:)=0;
for i=1:nz
    omega(ne,i)=(1/Re)*(Ex(ne,i)*((v(ne,i)-v(ne-1,i))/de)-Ey(ne,i)*((u(ne,i)-u(ne-1,i))/de));
    vort(ne,i)=Re*omega(ne,i);
   if i==1
       psi(1,ne,i)=U*Y(ne,i);
       psi(2,ne,i)=psi(1,ne,i);
   elseif i==nz
       psi(1,ne,i)=U*Y(ne,i);
       psi(2,ne,i)=psi(1,ne,i);
   else
       psi(1,ne,i)=U*Y(ne,i);
       psi(2,ne,i)=psi(1,ne,i);
   end
end
%Bottom
u(1,:)=0;
v(1,:)=0;
for i=1:nz
    omega(1,i)=(1/Re)*(Ex(1,i)*((v(2,i)-v(1,i))/de)-Ey(1,i)*((u(2,i)-u(1,i))/de));
    vort(1,i)=Re*omega(ne,i);
   if i==1
       psi(1,1,i)=U*Y(1,i);
       psi(2,1,i)=psi(1,1,i);
   elseif i==nz
       psi(1,1,i)=U*Y(1,i);
       psi(2,1,i)=psi(1,1,i);
   else
       psi(1,1,i)=U*Y(1,i);
       psi(2,1,i)=psi(1,1,i);
   end
end
%resultant velocity vector field
for j=1:ne
   for i=1:nz
      Vr(j,i)=((u(j,i)^2)+(v(j,i)^2))^0.5; 
   end
end
%plot
if menu==3
    X3=[X;X;X;X];
    Y3=[Y;(ly+0.01*de)*ones(size(Y))+Y;(2*ly+0.01*de)*ones(size(Y))+Y;(3*ly+0.01*de)*ones(size(Y))+Y];
    u3=[u;u;u;u];
    v3=[v;v;v;v];
    vr3=[Vr;Vr;Vr;Vr];
    figure(3)
    surf(X3,Y3,u3,'EdgeColor','none')
    view(0,90);
    colormap jet;
    title(['u at t=',num2str(T)]);
    xlabel('X');
    ylabel('Y');
    colorbar
    axis equal
    xlim([0 lx])
    ylim([dnorth 4*ly-dnorth])
    figure(4)
    surf(X3,Y3,v3,'EdgeColor','none')
    view(0,90);
    colormap jet;
    title(['v at t=',num2str(T)]);
    xlabel('X');
    ylabel('Y');
    colorbar
    axis equal
    xlim([0 lx])
    ylim([dnorth 4*ly-dnorth])
    figure(5)
    surf(X3,Y3,vr3,'EdgeColor','none')
    view(0,90);
    colormap jet;
    title(['Vr at t=',num2str(T)]);
    xlabel('X');
    ylabel('Y');
    colorbar
    axis equal
    xlim([0 lx])
    ylim([dnorth 4*ly-dnorth])
else
    figure(3)%resultant velocity (no flow)
    subplot(2,1,1)%Sub1
    surf(X,Y,Vr,'EdgeColor','none');
    view(0,90);
    colormap jet;
    title(['Vr at t=',num2str(T)]);
    xlabel('X');
    ylabel('Y');
    colorbar
    axis equal
    xlim([0 lx])
    ylim([0 ly])
    subplot(2,1,2)%Sub2
    quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),u(1:2:end,1:2:end),v(1:2:end,1:2:end))
    xlim([0,lx])
    ylim([0,ly])
    title(['Vector field at t=',num2str(T)])
    xlabel('X')
    ylabel('Y')
    axis equal
    xlim([0 lx])
    ylim([0 ly])
    %Vplot(1)=getframe(gcf);
%     figure(4)%vorticity & stream function
%     subplot(2,1,1)%Sub1
%     contour(X,Y,vort,70)
%     view(0,90);
%     colormap jet;
%     title(['\omega at t=',num2str(T)]);
%     xlabel('X');
%     ylabel('Y');
%     colorbar
%     pbaspect([lx,ly,1])
%     subplot(2,1,2)%Sub2
%     for i=1:nz
%        for j=1:ne
%            ppsi(j,i)=psi(1,j,i);
%        end
%     end
%     contour(X,Y,ppsi,70)
%     view(0,90);
%     colormap jet;
%     title(['\psi at t=',num2str(T)]);
%     xlabel('X');
%     ylabel('Y');
%     colorbar
%     pbaspect([lx,ly,1])
    % VSplot(1)=getframe(gcf);
    figure(5)%Axial and radial velocities
    subplot(2,1,1)%Sub1
    surf(X,Y,u,'EdgeColor','none')
    view(0,90);
    colormap jet;
    colorbar
    title(['u at t=',num2str(T)]);
    xlabel('X');
    ylabel('Y');
    axis equal
    xlim([0 lx])
    ylim([0 ly])
    subplot(2,1,2)%Sub2
    surf(X,Y,v,'EdgeColor','none')
    view(0,90);
    colormap jet;
    colorbar
    title(['v at t=',num2str(T)]);
    xlabel('X');
    ylabel('Y');
    axis equal
    xlim([0 lx])
    ylim([0 ly])
    % uvSurfPlots(1)=getframe(gcf);
end
 
%Main Loop
np=2;
count=0;
dPSI=1;
% dPSIp=1000;
if Re~=100
    if Re>100 && Re<=1000
        Rstp=5;
    elseif Re>1000 && Re<=10000
        Rstp=10;
    elseif Re>10000 && Re<=100000
        Rstp=50;
    end
    Rer=linspace(100,Re,Rstp);
    rc=1;
    Re=Rer(1,rc);
elseif Re==100
    Rstp=1;
    rc=Rstp+1;
    sc=9e-6;
end
while dPSI>sc
    if rc<Rstp
        if dPSI<6e-4
            rc=rc+1;
            Re=Rer(1,rc);
        end
    end
    
    %BC
    %Left
    u(:,1)=U;
    v(:,1)=0;
    for j=2:ne-1
       psi(1,j,1)=U*Y(j,1);
       psi(2,j,1)=psi(1,j,1);
       omega(j,1)=(1/Re)*(Zx(j,1)*((v(j,2)-v(j,1))/dz)-Zy(j,1)*((u(j,2)-u(j,1))/dz));
       vort(j,1)=Re*omega(j,1);
    end
    %Right
    u(:,nz)=u(:,nz-1);
    v(:,nz)=v(:,nz-1);
    for j=2:ne-1
       psi(1,j,nz)=psi(1,j,nz-1);
       psi(2,j,nz)=psi(1,j,nz);
       omega(j,nz)=(1/Re)*(Ex(j,nz)*((v(j+1,nz)-v(j-1,nz))/(2*de))-Ey(j,nz)*((u(j+1,nz)-u(j-1,nz))/(2*de)));
       vort(j,nz)=Re*omega(j,nz);
    end
    %Top
    u(ne,:)=0;
    v(ne,:)=0;
    for i=1:nz
        omega(ne,i)=(1/Re)*(Ex(ne,i)*((v(ne,i)-v(ne-1,i))/de)-Ey(ne,i)*((u(ne,i)-u(ne-1,i))/de));
        vort(ne,i)=Re*omega(ne,i);
       if i==1
           psi(1,ne,i)=U*Y(ne,i);
           psi(2,ne,i)=psi(1,ne,i);
       elseif i==nz
           psi(1,ne,i)=U*Y(ne,i);
           psi(2,ne,i)=psi(1,ne,i);
       else
           psi(1,ne,i)=U*Y(ne,i);
           psi(2,ne,i)=psi(1,ne,i);
       end
    end
    %Bottom
    u(1,:)=0;
    v(1,:)=0;
    for i=1:nz
        omega(1,i)=(1/Re)*(Ex(1,i)*((v(2,i)-v(1,i))/de)-Ey(1,i)*((u(2,i)-u(1,i))/de));
        vort(1,i)=Re*omega(ne,i);
       if i==1
           psi(1,1,i)=U*Y(1,i);
           psi(2,1,i)=psi(1,1,i);
       elseif i==nz
           psi(1,1,i)=U*Y(1,i);
           psi(2,1,i)=psi(1,1,i);
       else
           psi(1,1,i)=U*Y(1,i);
           psi(2,1,i)=psi(1,1,i);
       end
    end
    u(:,1)=U;
        %Omega
         for it=1:inner
           for j=2:ne-1
                for i=2:nz-1
                    A=(Zx(j,i)^2)+(Zy(j,i)^2);
                    B=(Ex(j,i)^2)+(Ey(j,i)^2);
                    C=2*Zx(j,i)*Ex(j,i)+2*Zy(j,i)*Ey(j,i);
                    D=((Zx(j,i+1)-Zx(j,i-1))/(X(j,i+1)-X(j,i-1)))+Zx(j,i)*((Zx(j,i+1)-Zx(j,i-1))/(2*dz))+Ex(j,i)*((Zx(j+1,i)-Zx(j-1,i))/(2*de))+((Zy(j+1,i)-Zy(j-1,i))/(Y(j+1,i)-Y(j-1,i)))+Zy(j,i)*((Zy(j,i+1)-Zy(j,i-1))/(2*dz))+Ey(j,i)*((Zy(j+1,i)-Zy(j-1,i))/(2*de));
                    E=((Ex(j,i+1)-Ex(j,i-1))/(X(j,i+1)-X(j,i-1)))+Ex(j,i)*((Ex(j,i+1)-Ex(j,i-1))/(2*de))+Zx(j,i)*((Ex(j,i+1)-Ex(j,i-1))/(2*dz))+((Ey(j+1,i)-Ey(j-1,i))/(Y(j+1,i)-Y(j-1,i)))+Ey(j,i)*((Ey(j+1,i)-Ey(j-1,i))/(2*de))+Zy(j,i)*((Ey(j,i+1)-Ey(j,i-1))/(2*dz));
                    omega(j,i)=(((A/(dz^2))*(omega(j,i+1)+omega(j,i-1)))+((B/(de^2))*(omega(j+1,i)+omega(j-1,i)))+((C/(4*dz*de))*(omega(j+1,i+1)-omega(j-1,i+1)-omega(j+1,i-1)+omega(j-1,i-1)))+((D/(2*dz))*(omega(j,i+1)-omega(j,i-1)))+((E/(2*de))*(omega(j+1,i)-omega(j-1,i)))-((1/(2*dz))*((u(j,i)*Zx(j,i))+(v(j,i)*Zy(j,i)))*(vort(j,i+1)-vort(j,i-1)))-((1/(2*de))*((u(j,i)*Ex(j,i))+(v(j,i)*Ey(j,i)))*(vort(j+1,i)-vort(j-1,i))))/((2*A/(dz^2))+(2*B/(de^2)));
                end
            end 
         end
        %Psi
        for j=2:ne-1
           for i=2:nz-1
              A=(Zx(j,i)^2)+(Zy(j,i)^2);
              B=(Ex(j,i)^2)+(Ey(j,i)^2);
              C=2*Zx(j,i)*Ex(j,i)+2*Zy(j,i)*Ey(j,i);
              D=((Zx(j,i+1)-Zx(j,i-1))/(X(j,i+1)-X(j,i-1)))+Zx(j,i)*((Zx(j,i+1)-Zx(j,i-1))/(2*dz))+Ex(j,i)*((Zx(j+1,i)-Zx(j-1,i))/(2*de))+((Zy(j+1,i)-Zy(j-1,i))/(Y(j+1,i)-Y(j-1,i)))+Zy(j,i)*((Zy(j,i+1)-Zy(j,i-1))/(2*dz))+Ey(j,i)*((Zy(j+1,i)-Zy(j-1,i))/(2*de));
              E=((Ex(j,i+1)-Ex(j,i-1))/(X(j,i+1)-X(j,i-1)))+Ex(j,i)*((Ex(j,i+1)-Ex(j,i-1))/(2*de))+Zx(j,i)*((Ex(j,i+1)-Ex(j,i-1))/(2*dz))+((Ey(j+1,i)-Ey(j-1,i))/(Y(j+1,i)-Y(j-1,i)))+Ey(j,i)*((Ey(j+1,i)-Ey(j-1,i))/(2*de))+Zy(j,i)*((Ey(j,i+1)-Ey(j,i-1))/(2*dz));
              psi(2,j,i)=((dt/Re)*(((A/(dz^2))*(psi(1,j,i+1)-2*psi(1,j,i)+psi(1,j,i-1)))+((B/(de^2))*(psi(1,j+1,i)-2*psi(1,j,i)+psi(1,j-1,i)))+((C/(4*dz*de))*(psi(1,j+1,i+1)-psi(1,j-1,i+1)-psi(1,j+1,i-1)+psi(1,j-1,i-1)))+((D/(2*dz))*(psi(1,j,i+1)-psi(1,j,i-1)))+((E/(2*de))*(psi(1,j+1,i)-psi(1,j-1,i)))))+(dt*omega(j,i))+psi(1,j,i);
           end
        end
        dPSI=max(max(abs((psi(2,:,:)-psi(1,:,:))./psi(2,:,:))));%max change of stream function at each time step
%         if (dPSIp/dPSI)>=0.9
%             dPSIp=dPSI;
%         elseif (dPSIp/dPSI)<0.9
%             dPSI=1e-7;
%         end
        psi(1,:,:)=psi(2,:,:);
        %Vorticity
        for j=2:ne-1
           for i=2:nz-1
              A=(Zx(j,i)^2)+(Zy(j,i)^2);
              B=(Ex(j,i)^2)+(Ey(j,i)^2);
              C=2*Zx(j,i)*Ex(j,i)+2*Zy(j,i)*Ey(j,i);
              D=((Zx(j,i+1)-Zx(j,i-1))/(X(j,i+1)-X(j,i-1)))+Zx(j,i)*((Zx(j,i+1)-Zx(j,i-1))/(2*dz))+Ex(j,i)*((Zx(j+1,i)-Zx(j-1,i))/(2*de))+((Zy(j+1,i)-Zy(j-1,i))/(Y(j+1,i)-Y(j-1,i)))+Zy(j,i)*((Zy(j,i+1)-Zy(j,i-1))/(2*dz))+Ey(j,i)*((Zy(j+1,i)-Zy(j-1,i))/(2*de));
              E=((Ex(j,i+1)-Ex(j,i-1))/(X(j,i+1)-X(j,i-1)))+Ex(j,i)*((Ex(j,i+1)-Ex(j,i-1))/(2*de))+Zx(j,i)*((Ex(j,i+1)-Ex(j,i-1))/(2*dz))+((Ey(j+1,i)-Ey(j-1,i))/(Y(j+1,i)-Y(j-1,i)))+Ey(j,i)*((Ey(j+1,i)-Ey(j-1,i))/(2*de))+Zy(j,i)*((Ey(j,i+1)-Ey(j,i-1))/(2*dz));
              vort(j,i)=(((-1*A/(dz^2))*(psi(1,j,i+1)-2*psi(1,j,i)+psi(1,j,i-1)))-((B/(de^2))*(psi(1,j+1,i)-2*psi(1,j,i)+psi(1,j-1,i)))-((C/(4*dz*de))*(psi(1,j+1,i+1)-psi(1,j-1,i+1)-psi(1,j+1,i-1)+psi(1,j-1,i-1)))-((D/(2*dz))*(psi(1,j,i+1)-psi(1,j,i-1)))-((E/(2*de))*(psi(1,j+1,i)-psi(1,j-1,i))));
           end
        end
        %Velocities
        for j=2:ne-1
           for i=2:nz-1
              u(j,i)=((Zy(j,i)/(2*dz))*(psi(1,j,i+1)-psi(1,j,i-1)))+((Ey(j,i)/(2*de))*(psi(1,j+1,i)-psi(1,j-1,i)));
              v(j,i)=((Zx(j,i)/(-2*dz))*(psi(1,j,i+1)-psi(1,j,i-1)))-((Ex(j,i)/(2*de))*(psi(1,j+1,i)-psi(1,j-1,i)));
           end
        end
        count=count+1;
        T=T+dt;
 
        %Plot
        if count==q
            %resultant velocity vector field
            for j=1:ne
                for i=1:nz
                    Vr(j,i)=((u(j,i)^2)+(v(j,i)^2))^0.5; 
                end
            end
            if menu==3
                X3=[X;X;X;X];
                Y3=[Y;(ly+0.01*de)*ones(size(Y))+Y;(2*ly+0.01*de)*ones(size(Y))+Y;(3*ly+0.01*de)*ones(size(Y))+Y];
                u3=[u;u;u;u];
                v3=[v;v;v;v];
                vr3=[Vr;Vr;Vr;Vr];
                figure(3)
                surf(X3,Y3,u3,'EdgeColor','none')
                view(0,90);
                colormap jet;
                title(['u at t=',num2str(T)]);
                xlabel('X');
                ylabel('Y');
                colorbar
                axis equal
                xlim([0 lx])
                ylim([dnorth 4*ly-dnorth])
                figure(4)
                surf(X3,Y3,v3,'EdgeColor','none')
                view(0,90);
                colormap jet;
                title(['v at t=',num2str(T)]);
                xlabel('X');
                ylabel('Y');
                colorbar
                axis equal
                xlim([0 lx])
                ylim([dnorth 4*ly-dnorth])
                figure(5)
                surf(X3,Y3,vr3,'EdgeColor','none')
                view(0,90);
                colormap jet;
                title(['Vr at t=',num2str(T)]);
                xlabel('X');
                ylabel('Y');
                colorbar
                axis equal
                xlim([0 lx])
                ylim([dnorth 4*ly-dnorth])
            else
                figure(3)%resultant velocity (no flow)
                subplot(2,1,1)%Sub1
                surf(X,Y,Vr,'EdgeColor','none');
                view(0,90);
                colormap jet;
                title(['Vr at t=',num2str(T)]);
                xlabel('X');
                ylabel('Y');
                colorbar
                axis equal
                xlim([0 lx])
                ylim([0 ly])
                subplot(2,1,2)%Sub2
                quiver(X(1:2:end,1:4:end),Y(1:2:end,1:4:end),u(1:2:end,1:4:end),v(1:2:end,1:4:end))
                xlim([0,lx])
                ylim([0,ly])
                title(['Vector field at t=',num2str(T)])
                xlabel('X')
                ylabel('Y')
                axis equal
                xlim([0 lx])
                ylim([0 ly])
                %Vplot(1)=getframe(gcf);
%                 figure(4)%vorticity & stream function
%                 subplot(2,1,1)%Sub1
%                 contour(X,Y,vort,70)
%                 view(0,90);
%                 colormap jet;
%                 title(['\omega at t=',num2str(T)]);
%                 xlabel('X');
%                 ylabel('Y');
%                 colorbar
%                 pbaspect([lx,ly,1])
%                 subplot(2,1,2)%Sub2
%                 for i=1:nz
%                    for j=1:ne
%                        ppsi(j,i)=psi(1,j,i);
%                    end
%                 end
%                 contour(X,Y,ppsi,70)
%                 view(0,90);
%                 colormap jet;
%                 title(['\psi at t=',num2str(T)]);
%                 xlabel('X');
%                 ylabel('Y');
%                 colorbar
%                 pbaspect([lx,ly,1])
                % VSplot(1)=getframe(gcf);
                figure(5)%Axial and radial velocities
                subplot(2,1,1)%Sub1
                surf(X,Y,u,'EdgeColor','none')
                view(0,90);
                colormap jet;
                colorbar
                title(['u at t=',num2str(T)]);
                xlabel('X');
                ylabel('Y');
                axis equal
                xlim([0 lx])
                ylim([0 ly])
                subplot(2,1,2)%Sub2
                surf(X,Y,v,'EdgeColor','none')
                view(0,90);
                colormap jet;
                colorbar
                title(['v at t=',num2str(T)]);
                xlabel('X');
                ylabel('Y');
                axis equal
                xlim([0 lx])
                ylim([0 ly])
                % uvSurfPlots(1)=getframe(gcf);
            end
            count=0;
            np=np+1;
        end
end
 
if menu==1
   ind=[4,12,81,161,169];
   figure(6)
   for k=1:5
       subplot(1,5,k)
       plot(u(:,ind(1,k)),Y(:,ind(1,k)),'LineWidth',2)
       title(['u at x=',num2str(X(1,ind(1,k)))])
       xlabel('u')
       ylabel('Y')
   end
   set(gcf, 'Position',  [10, 10, 1500, 300])
   
elseif menu==2
   ind=[41,57,73,89,105,121];
   figure(6)
   for k=1:6
       subplot(1,6,k)
       plot(u(:,ind(1,k)),Y(:,ind(1,k)),'LineWidth',2)
       title(['u at x=',num2str(X(1,ind(1,k)))])
       xlabel('u')
       ylabel('Y')
   end
   set(gcf, 'Position',  [10, 10, 1800, 300])
   
elseif menu==3
    Cp=ones(size(u))-u.^2;
    Cl=trapz(X(3,3:gd),(Cp(py-2,3:gd)-Cp(3,3:gd)));
    figure(6)
%     plot(X(3,3:gd),Cp(py-2,3:gd)-Cp(3,3:gd),'LineWidth',2)
    area(X(3,3:gd),Cp(py-2,3:gd)-Cp(3,3:gd))
    title(['Cl=',num2str(Cl)])
    xlabel('X')
    ylabel('Cpl-Cpu')
    
end
    elseif op==3
        op=NaN;
%         disp('Option3')
        clear variables
%Input parameters
nx=input('Grid density (points/unit length): ');%grid points in x-direction
ny=nx;%grid points in y-direction (same as x, dx=dy)
U=input('Lid velocity oscillation amplitude (U): ');%Top wall amplitude velocity
an=input('Angular velocity (w=n*pi), n: ');%angular velocity
Re=input('Reynolds number: ');%Reynolds number
cycles=input('Number of cycles: ');%number of cycles
ts=input(['Total time is: ',num2str(cycles*2/an),'s, Specify timesteps (','recommended: ',num2str((cycles*2/an)/(0.5*Re*((1/nx)^1)*((1/ny)^1))),'<<): ']);%number of time steps
q=input(['ts/100=',num2str(ts/100),' Plot every: ']);

w=an*pi;%angular velocity
X=linspace(0,1,nx);%x-coordinate points
Y=X;%y-coordinate points
dx=X(1,2)-X(1,1);%dx
dy=dx;%dy
T=linspace(0,cycles*2/an,ts);%time coordinate points
uc=zeros(1,length(T));vc=uc;%center point velocities
dt=T(1,2)-T(1,1);%dt
u=zeros(ny,nx);%u-velocity matrix
v=u;%v-velocity matrix
vort=u;%vorticity matrix
psi=zeros(2,ny,nx);%stream matrix
omega=u;%Omega matrix
Vr=u;%resultant velocity matrix

%BC's
u(1,:)=0; v(1,:)=0; psi(1,1,:)=0;
v(ny,:)=0; psi(1,ny,:)=0;
u(:,1)=0; v(:,1)=0; psi(1,:,1)=0;
u(:,nx)=0; v(:,nx)=0; psi(1,:,nx)=0;
for i=1:nx
   omega(1,i)=(-1/Re)*((u(2,i)-u(1,i))/dy);
   vort(1,i)=Re*omega(1,i);
   u(ny,i)=U*cos(w*T(1,1));
   omega(ny,i)=(-1/Re)*((u(ny,i)-u(ny-1,i))/dy);
   vort(ny,i)=Re*omega(ny,i);
end
for j=1:ny
    omega(j,1)=(1/Re)*((v(j,2)-v(j,1))/dx);
    vort(j,1)=Re*omega(j,1);
    omega(j,nx)=(1/Re)*((v(j,nx)-v(j,nx-1))/dx);
    vort(j,nx)=Re*omega(j,nx);
end
up=u;
vp=v;

menui=2;
while menui~=0 && menui~=1
    menui=input('Start from steady state solution? \n [1]: Yes \n [0]: No \n Selection: ');
end
if menui==1
    urms=1;vrms=1;it=1;
    while (urms>1e-6 || vrms>1e-6) && it<1000
        for j=2:ny-1
           for i=2:nx-1
              omega(j,i)=(((omega(j,i+1)+omega(j,i-1))/(dx^2))+((omega(j+1,i)+omega(j-1,i))/(dy^2))-(0.5*u(j,i)*(vort(j,i+1)-vort(j,i-1))/dx)-(0.5*v(j,i)*(vort(j+1,i)-vort(j-1,i))/dy))/((2/dx^2)+(2/dy^2));
           end
        end
        
        for j=2:ny-1
           for i=2:nx-1
               psi(1,j,i)=(((psi(1,j,i+1)+psi(1,j,i-1))/(dx^2))+((psi(1,j+1,i)+psi(1,j-1,i))/(dy^2))+(Re*omega(j,i)))/((2/(dx^2))+(2/(dy^2)));
           end
        end
        
        for j=2:ny-1
           for i=2:nx-1
              vort(j,i)=-((psi(1,j,i+1)-2*psi(1,j,i)+psi(1,j,i-1))/(dx^2))-((psi(1,j+1,i)-2*psi(1,j,i)+psi(1,j-1,i))/(dy^2));
           end
        end
        
        sumu=0; sumv=0;
        for j=2:ny-1
           for i=2:nx-1
              u(j,i)=((psi(1,j+1,i)-psi(1,j-1,i))/(2*dy)); 
              sumu=sumu+((u(j,i)-up(j,i))^2);
              v(j,i)=-((psi(1,j,i+1)-psi(1,j,i-1))/(2*dx));
              sumv=sumv+((v(j,i)-vp(j,i))^2);
           end
        end
        urms=(sumu/((nx-1)*(ny-1)))^0.5;
        vrms=(sumv/((nx-1)*(ny-1)))^0.5;
        up=u;
        vp=v;
        
        for i=1:nx
           omega(1,i)=(-1/Re)*((u(2,i)-u(1,i))/dy);
           vort(1,i)=Re*omega(1,i);
           omega(ny,i)=(-1/Re)*((u(ny,i)-u(ny-1,i))/dy);
           vort(ny,i)=Re*omega(ny,i);
        end
        for j=1:ny
            omega(j,1)=(1/Re)*((v(j,2)-v(j,1))/dx);
            vort(j,1)=Re*omega(j,1);
            omega(j,nx)=(1/Re)*((v(j,nx)-v(j,nx-1))/dx);
            vort(j,nx)=Re*omega(j,nx);
        end
        it=it+1;
    end
end

for j=1:ny
   for i=1:nx
      Vr(j,i)=((u(j,i)^2)+(v(j,i)^2))^0.5; 
   end
end

figure(1)
subplot(2,2,1)
surfc(X,Y,Vr,'EdgeColor','none');
view(0,90);
colormap jet;
title(['Vr at t=',num2str(T(1,1))]);
xlabel('X');
ylabel('Y');
colorbar
pbaspect([1,1,1])
subplot(2,2,2)
colormap jet
contour(X,Y,Vr)
hold on
streamslice(X,Y,u,v)
% xlim([0,1])
% ylim([0,1])
%quiver(X,Y,u,v,'AutoScale','on', 'AutoScaleFactor', 4)
hold off
title(['Velocity vector field at t=',num2str(T(1,1))])
xlabel('X')
ylabel('Y')
colorbar
pbaspect([1,1,1])
subplot(2,2,3)
plot(u(:,(nx+1)/2),Y)
title('Vertical center line u-velocity')
xlabel('u-velocity')
ylabel('Y')
pbaspect([1,1,1])
subplot(2,2,4)
plot(X,v((ny+1)/2,:))
title('Horizontal center line v-velocity')
xlabel('X')
ylabel('v-velocity')
pbaspect([1,1,1])
Vplot(1)=getframe(gcf);

uc(1,1)=u((ny+1)/2,(nx+1)/2);
vc(1,1)=v((ny+1)/2,(nx+1)/2);

np=2;
%count=1;
for k=2:ts
    %BC
    for i=1:nx
       omega(1,i)=(-1/Re)*((u(2,i)-u(1,i))/dy);
       vort(1,i)=Re*omega(1,i);
       u(ny,i)=U*cos(w*T(1,k));
       omega(ny,i)=(-1/Re)*((u(ny,i)-u(ny-1,i))/dy);
       vort(ny,i)=Re*omega(ny,i);
    end
    for j=1:ny
        omega(j,1)=(1/Re)*((v(j,2)-v(j,1))/dx);
        vort(j,1)=Re*omega(j,1);
        omega(j,nx)=(1/Re)*((v(j,nx)-v(j,nx-1))/dx);
        vort(j,nx)=Re*omega(j,nx);
    end
    %Omega
    for j=2:ny-1
       for i=2:nx-1
          omega(j,i)=(((omega(j,i+1)+omega(j,i-1))/(dx^2))+((omega(j+1,i)+omega(j-1,i))/(dy^2))-(0.5*u(j,i)*(vort(j,i+1)-vort(j,i-1))/dx)-(0.5*v(j,i)*(vort(j+1,i)-vort(j-1,i))/dy))/((2/dx^2)+(2/dy^2));
       end
    end
    %psi
    for j=2:ny-1
       for i=2:nx-1
          psi(2,j,i)=(dt/Re)*(((psi(1,j,i+1)-2*psi(1,j,i)+psi(1,j,i-1))/(dx^2))+((psi(1,j+1,i)-2*psi(1,j,i)+psi(1,j-1,i))/(dy^2)))+dt*omega(j,i)+psi(1,j,i);
       end
    end
    psi(1,:,:)=psi(2,:,:);
    %vorticity
    for j=2:ny-1
       for i=2:nx-1
          vort(j,i)=-((psi(1,j,i+1)-2*psi(1,j,i)+psi(1,j,i-1))/(dx^2))-((psi(1,j+1,i)-2*psi(1,j,i)+psi(1,j-1,i))/(dy^2));
       end
    end
    %velocities
    for j=2:ny-1
       for i=2:nx-1
          u(j,i)=((psi(1,j+1,i)-psi(1,j-1,i))/(2*dy)); 
          v(j,i)=-((psi(1,j,i+1)-psi(1,j,i-1))/(2*dx));
       end
    end
    count=count+1;
    
    uc(1,k)=u((ny+1)/2,(nx+1)/2);
    vc(1,k)=v((ny+1)/2,(nx+1)/2);
    
    if count==q
       %Resultant velocity and plot
       for j=1:ny
          for i=1:nx
             Vr(j,i)=((u(j,i)^2)+(v(j,i)^2))^0.5; 
          end
       end
       %figure(1)
       subplot(2,2,1)
       surfc(X,Y,Vr,'EdgeColor','none');
       view(0,90);
       colormap jet;
       title(['Vr at t=',num2str(T(1,k))]);
       xlabel('X');
       ylabel('Y');
       colorbar
       pbaspect([1,1,1])
       subplot(2,2,2)
       colormap jet
       contour(X,Y,Vr)
       hold on
       streamslice(X,Y,u,v)
%        xlim([0,1])
%        ylim([0,1])
%        quiver(X,Y,u,v,'AutoScale','on', 'AutoScaleFactor', 4)
       hold off
       title(['Velocity vector field at t=',num2str(T(1,k))])
       xlabel('X')
       ylabel('Y')
       colorbar
       pbaspect([1,1,1])
       subplot(2,2,3)
       plot(u(:,(nx+1)/2),Y)
       title('Vertical center line u-velocity')
       xlabel('u-velocity')
       ylabel('Y')
       pbaspect([1,1,1])
       subplot(2,2,4)
       plot(X,v((ny+1)/2,:))
       title('Horizontal center line v-velocity')
       xlabel('X')
       ylabel('v-velocity')
       pbaspect([1,1,1])
       Vplot(np)=getframe(gcf);
       count=1;
       np=np+1;
    end
    
end

menuV=2;
while menuV~=0 && menuV~=1
    menuV=input('Make Video? \n [1]: Yes \n [0]: No \n Selection: ');
end
if menuV==1
    fps=input(['Available frames: ',num2str(length(Vplot)),', Specify fps:']);
    videoVr=VideoWriter('Vr.avi');
    videoVr.FrameRate=fps;
    videoVr.Quality=100;
    open(videoVr)
    writeVideo(videoVr,Vplot)
    close(videoVr)
end

figure(2)
yyaxis left
plot(T,uc)
yyaxis right
plot(T,vc)
yyaxis left
title('Center point velocity in time')
xlabel('Time')
ylabel('u-velocity')
yyaxis right
ylabel('v-velocity')
    elseif op==4
        op=NaN;
%         disp('Option4')
        clear variables

%MAKE SURE THIS PARAMETERS ARE CORRECT FOR THE VECTOR FIELD DAT FILES!
delimiterIn = ' ';%columns delimiter (one space)
headerlinesIn = 3;%3 headers before data
ni=1;%Number of images = 1 (no stitching)

nf=input('How many files per image?: ');%Number of vector fields available
cfold=uigetdir('C:\Users\rodrigo\Documents\','Select Folder Containing Vector Fields')
mainfolder=cd(cfold);
for g=1:nf
    i_strPadded = sprintf( '%05d', g );%create file number
    filename=['B',i_strPadded,'.dat'];%create file name
    A = importdata(filename,delimiterIn,headerlinesIn);%Import data to A (structure variable)

    if m==1 && g==1%Only the first time a file is read, matrix sizes are calculated and variables are created

       diff=0;c=1;n=0;
       while diff==0
          n=n+1;
          diff=A.data(c,2)-A.data(c+1,2);
          c=c+1;
       end
       nx=n;%matrix size (nx) of each image
       ny=size(A.data,1)/nx;%matrix size (ny) of each image

       xi=zeros(nx,1);%x-coordinates (assuming the PIV processing is the same for all images)
       yi=zeros(ny,1);%y-coordinates 
       Avx=zeros(ny,nx,nf,ni);%A-matrix for axial velocity component (Vx)
       Avy=Avx;%A-matrix for radial velocity component (Vy)
       Avr=Avx;%A-matrix for resultant velocity component (Vr)

       for z=1:ny
           yi(z,1)=A.data((ny*nx)-nx*(z-1),2);
       end
       for z=1:nx
          xi(z,1)=A.data(z,1);
       end
       dx=xi(2,1)-xi(1,1);%dx increment in x-coordinates
       dy=yi(2,1)-yi(1,1);%dy increment in y-coordinates

    end

    %Filling up matrices
    r=1;
    for j=1:ny
       for i=1:nx
          Avx((ny+1)-j,i,g,m)=A.data(r,3);%filling up Vx
          Avy((ny+1)-j,i,g,m)=A.data(r,4);%filling up Vy
          Avr((ny+1)-j,i,g,m)=(((Avx((ny+1)-j,i,g,m))^2)+((Avy((ny+1)-j,i,g,m))^2))^0.5;%filling up Vr
          r=r+1; 
       end
    end

end
cd (mainfolder);%Go back to main folder

Vx=Avx;
Vy=Avy;

ok=2;it=1;
while ok~=0
    disp(' ')
    disp(['Iteration: ', num2str(it)])
    xl=input('Number of columns at left end to be removed: ');
    xr=input('Number of columns at right end to be removed: ');
    yt=input('Number of top rows to be removed: ');
    yb=input('Number of bottom rows to be removed: ');

    aVx=Vx(yb+1:ny-yt,xl+1:nx-xr,:);
    aVy=Vy(yb+1:ny-yt,xl+1:nx-xr,:);
    Y=yi(yb+1:ny-yt,1);
    anx=size(aVx,2);
    X=zeros(1,anx);
    for i=1:anx-1
        X(1,i+1)=i*dx;
    end
    figure(it)
    subplot(2,1,1)
    surf(X,Y,mean(aVx,3),'EdgeColor','none')
    axis equal
    colormap jet
    view(0,90)
    title('Vx')
    colorbar
    subplot(2,1,2)
    surf(X,Y,mean(aVy,3),'EdgeColor','none')
    axis equal
    colormap jet
    view(0,90)
    title('Vy')
    colorbar
    set(gcf, 'Position',  [10, 10, 1200, 1000])

    ok=input('Edit more? [1]:YES [0]:NO \n Selection: ');
    it=it+1;
end
Vx=aVx;
Vy=aVy;
nx=anx;
ny=size(Vx,1);

fname=input('Input name: ','s');
save([fname,'.mat'],'Vx','Vy','X','Y','nx','ny')
clearvars -except fname
load([fname,'.mat'])

%POD
nf=size(Vx,3);
Avx=zeros(nf,nx*ny);
Avy=Avx;
Avx_mean_sub=Avx;
Avy_mean_sub=Avx;

for k=1:nf
    r=1;
    for j=1:ny
       for i=1:nx
           Avx(k,r)=Vx(j,i,k);
           Avy(k,r)=Vy(j,i,k);
           r=r+1;
       end
    end
end

Avx_col_mean=mean(Avx,1);
Avy_col_mean=mean(Avy,1);

for i=1:nx*ny
   for j=1:nf
      Avx_mean_sub(j,i)=Avx(j,i)-Avx_col_mean(1,i);
      Avy_mean_sub(j,i)=Avy(j,i)-Avy_col_mean(1,i);
   end
end

[Uvx,Svx,Vvx]=svd(Avx_mean_sub,'econ');
[Uvy,Svy,Vvy]=svd(Avy_mean_sub,'econ');

lambdax=zeros(nf,1);
lambday=lambdax;
for i=1:nf
   lambdax(i,1)=Svx(i,i)^2;
   lambday(i,1)=Svy(i,i)^2;
end

PECx=(100/sum(lambdax))*lambdax;
PECy=(100/sum(lambday))*lambday;

d=0;rx=0;
while d<0.9
    rx=rx+1;
    d=sum(lambdax(1:rx,1))/sum(lambdax);
end

d=0;ry=0;
while d<0.9
    ry=ry+1;
    d=sum(lambday(1:ry,1))/sum(lambday);
end

PODmodeVx=zeros(ny,nx);
PODmodeVy=PODmodeVx;

for m=1:20
   r=1;
   for j=1:ny
      for i=1:nx
         PODmodeVx(j,i)=Vvx(r,m);
         PODmodeVy(j,i)=Vvy(r,m);
         r=r+1;
      end
   end
   
   figure(1)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVx,'EdgeColor','none')
   axis equal
   xlabel('X/De')
   ylabel('Y/De')
   colormap jet
   view(0,90)
   colorbar
   title(['Vx POD Mode ',num2str(m)])
   set(gcf, 'Position',  [10, 10, 1400, 400])
   Modex(m)=getframe(gcf);
   
   figure(2)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVy,'EdgeColor','none')
   axis equal
   xlabel('X/De')
   ylabel('Y/De')
   colormap jet
   view(0,90)
   colorbar
   title(['Vy POD Mode ',num2str(m)])
   set(gcf, 'Position',  [10, 10, 1400, 400])
   Modey(m)=getframe(gcf);
   
end

% fps=1;
% videoMx=VideoWriter([fname,'VxPOD.avi']);
% videoMx.FrameRate=fps;
% videoMx.Quality=100;
% open(videoMx)
% writeVideo(videoMx,Modex)
% close(videoMx)
% videoMy=VideoWriter([fname,'VyPOD.avi']);
% videoMy.FrameRate=fps;
% videoMy.Quality=100;
% open(videoMy)
% writeVideo(videoMy,Modey)
% close(videoMy)

figure(1)
for h=1:4
   r=1;
   for j=1:ny
      for i=1:nx
          PODmodeVx(j,i)=Vvx(r,h);
          PODmodeVy(j,i)=Vvy(r,h);
          r=r+1;
      end
   end
   subplot(4,2,(h*2)-1)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVx,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vx POD Mode ',num2str(h)])
   subplot(4,2,h*2)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVy,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vy POD Mode ',num2str(h)])
   set(gcf, 'Position',  [10, 10, 1500, 1000])
end

figure(2)
m=h+1;
for h=1:4
   r=1;
   for j=1:ny
      for i=1:nx
          PODmodeVx(j,i)=Vvx(r,m);
          PODmodeVy(j,i)=Vvy(r,m);
          r=r+1;
      end
   end
   subplot(4,2,(h*2)-1)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVx,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vx POD Mode ',num2str(m)])
   subplot(4,2,h*2)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVy,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vy POD Mode ',num2str(m)])
   set(gcf, 'Position',  [10, 10, 1500, 1000])
   m=m+1;
end

figure(3)
for h=1:4
   r=1;
   for j=1:ny
      for i=1:nx
          PODmodeVx(j,i)=Vvx(r,m);
          PODmodeVy(j,i)=Vvy(r,m);
          r=r+1;
      end
   end
   subplot(4,2,(h*2)-1)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVx,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vx POD Mode ',num2str(m)])
   subplot(4,2,h*2)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVy,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vy POD Mode ',num2str(m)])
   set(gcf, 'Position',  [10, 10, 1500, 1000])
   m=m+1;
end

figure(4)
for h=1:4
   r=1;
   for j=1:ny
      for i=1:nx
          PODmodeVx(j,i)=Vvx(r,m);
          PODmodeVy(j,i)=Vvy(r,m);
          r=r+1;
      end
   end
   subplot(4,2,(h*2)-1)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVx,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vx POD Mode ',num2str(m)])
   subplot(4,2,h*2)
   surf((1/80.416)*X,(1/80.416)*Y,PODmodeVy,'EdgeColor','none')
   axis equal
   colormap jet
   view(0,90)
   title(['Vy POD Mode ',num2str(m)])
   set(gcf, 'Position',  [10, 10, 1500, 1000])
   m=m+1;
end

figure (5)
subplot(1,2,1)
plot(PECx)
xlabel('POD Mode')
ylabel('Percent Energy Content')
title('Vx POD Energy Content')
subplot(1,2,2)
plot(PECy)
xlabel('POD Mode')
ylabel('Percent Energy Content')
title('Vy POD Energy Content')
set(gcf, 'Position',  [10, 10, 900, 400])

    elseif op==0
        disp('Goodbye')
    end
end