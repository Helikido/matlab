% Final exam code 

clear all;
clc; 
close all; 

% User variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Mesh size 
Nx=50;
Ny=50;

% Flow properties 
Re=100;
U=10;

% Boundary conditions 
u_left=0.0; 
u_right=0.0;
u_top=U;
u_bottom=0.0;

v_left=0.0;
v_right=0.0;
v_top=0.0;
v_bottom=0.0;

x_min=-0.5;
x_max=0.5;
y_min=-0.5;
y_max=0.5;

% Non user variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Lx=x_max-x_min;
Ly=y_max-y_min;
N=Nx*Ny;
N_u =(Nx-1)*Ny;
N_v=Nx*(Ny-1);
dx=Lx/Nx;
dy=Ly/Ny;
Nu=U/Re;

% Descritize x & y domains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x(1)=x_min+dx/2;
for i=2:Nx
    x(i)=x(i-1)+dx;
end

y(1)=y_min+dy/2;
for j=2:Ny
    y(j)=y(j-1)+dy;
end 

% Create mesh grid ~~~~~~~~~~~~~~~~~~~~~~~~
[X,Y]=meshgrid(x,y);

% Initialize guess values for p,u,v as zeros ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
u_guess=zeros(N_u,1);
v_guess=zeros(N_v,1);
p_guess=zeros(N,1);

% Initialize n constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~
residual=1.0;
n=0;

% Initialize main while loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while((residual>=1E-3) && (n<=10))
    
    % Initialize coefficients for x-momentum equation
    ae=zeros(N_u,1);
    aw=zeros(N_u,1);
    as=zeros(N_u,1);
    an=zeros(N_u,1);
    su=zeros(N_u,1);
    sp=zeros(N_u,1);

    % X-momentum equation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for j=1:Ny % Build in terms of the xmomentum equation u nodes only
        for i=1:Nx-1 % Offset by one for u in the xmomentum equation
            cell=(j-1)*(Nx-1)+i; % For the u cell cv 
            
            % Dealing with Fe, Fw, Fn, Fs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % For Fe nodes not at boundary
            if (i<Nx-1) 
                fe=0.5*(u_guess(cell)+u_guess(cell+1))*dy;
            end
            if (i==Nx-1) % Don't blindly take average of the last boundary column
                fe=0.5*(u_guess(cell)+u_right)*dy;
            end 
            
            % For Fw nodes not at boundary 
            if (i>1)
                fw=0.5*(u_guess(cell)+u_guess(cell-1))*dy;
            end
            if (i==1) % i=1 dont blindly take average of last column
                fw=0.5*(u_guess(cell)+u_left)*dy;
            end

            cell_v=cell+j-1; % For the v cell cv
            
            if (j<Ny) % for all the cells except for the top most row because it is given by the bioundary condition
                fn=0.5*(v_guess(cell_v)+v_guess(cell_v+1))*dx; % for all the cells 
            end
            if (j==Ny)
                fn=v_top*dx; % for the top row 
            end 
        
            if (j>1) % for all cells except for the bottom most row because its given by the boundary condition
                fs=0.5*(v_guess(cell_v-Nx)+v_guess(cell_v+1-Nx))*dx;
            end
            if (j==1)
                fs=v_bottom*dx; % bottom boundary condition 
            end
            
            % initializing diffusion terms and solving them ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            de=Nu/dx*dy;
            dw=Nu/dx*dy;
            
            if (j<Ny) % Same thing for north 
                dn=Nu/dy*dx;
            end
            if (j==Ny)
                dn=2*Nu/dy*dx; % only for the top boundary cells, because of the distance between the node and boundary is DX/2
            end
            
            if (j>1) % Same thing for u_top
                ds=Nu/dy*dx;
            end
            if (j==1)
                ds=2*Nu/dy*dx; %o nly for the bottom boundary cells, because of the distance between the node and boundary is DX/2
            end
            
            % Sets up all the coefficients, ap is not because its done already, its jsu_top p the sume of all coefficients and apply UDS logic
            aw(cell)=dw+max(fw,0);
            ae(cell)=de+max(0,-fe);
            an(cell)=dn+max(0,-fn);
            as(cell)=ds+max(fs,0);
            
            % Dealing with the dp/dx term ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            cell_p=cell+j-1; % For the P cells, the reason why it's not as usual is because of the Nx in the i loop.
            
            su(cell)=(p_guess(cell_p)-p_guess(cell_p+1))/dx*dy*dx; % pressure gradient is a source term, it goes to b ~~~~~~~~~~~~~~~~~~~~~~~~~
            sp(cell)=0.0; % internal cells ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            if (j==1) % Bottom boundary 
                su(cell)=su(cell)+as(cell)*u_bottom ;
            end
            
            if (j==Ny) % Top boundary 
                su(cell)=su(cell)+an(cell)*u_top;
            end
            
            if (i==1) % Left boundary 
                su(cell)=su(cell)+aw(cell)*u_left;
            end
            
            if (i==Nx) % Right boundary  
                su(cell)=su(cell)+ae(cell)*u_right;
            end
            
            ap(cell)=aw(cell)+ae(cell)+as(cell)+an(cell)-sp(cell);
            
            du(cell)=1.0/dx/ap(cell); % d coefficient to be used in pressure correction equation
        end
    end
        % Setting up the matrix

    A = zeros(N_u,N_u);
    for cell=1:N_u
            A(cell,cell)=ap(cell);
        if ((cell-1)>0)
            A(cell,cell-1)=-aw(cell);
        end
        if ((cell+1)<=N_u)
            A(cell,cell+1)=-ae(cell);
        end
        if ((cell-Nx+1)>0)
            A(cell,cell-Nx+1)=-as(cell);
        end
        if ((cell+Nx-1)<=N_u)
            A(cell,cell+Nx-1)=-an(cell);
        end
    end

    % Setting up the RHS
    b=zeros(N_u,1);
    for cell=1:N_u
        b(cell)=su(cell);
    end

    % Solve for u_guess
    u_guess_updated = A\b ;

    % Setting up coefficients for discretized y-momentum equation
    ae = zeros(N_v,1);
    aw = zeros(N_v,1);
    as = zeros(N_v,1);
    an = zeros(N_v,1);
    su = zeros(N_v,1);
    sp = zeros(N_v,1);
    
    % Y-Momentum equation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for j=1:Ny-1
        for i=1:Nx
            cell=(j-1)*(Nx)+i;
            
            cell_u=cell-j+1;
            if (i<Nx)
                fe=0.5*(u_guess(cell_u)+u_guess(cell_u+Nx-1))*dy;
            end 
            if (i==Nx)
                fe=u_right*dy;
            end
            
            if (i>1)
                fw =0.5*(u_guess(cell_u-1)+u_guess(cell_u-1+Nx-1))*dy;
            end
            if (i==1)
                fw=u_left*dy;
            end
            
            
            if (j<Ny-1)
                fn=0.5*(v_guess(cell)+v_guess(cell+Nx))*dx;
            end
            if (j==(Ny-1))
                fn=0.5*(v_guess(cell)+v_top)*dx;
            end
            
            if (j>1)
                fs=0.5*(v_guess(cell)+v_guess(cell-Nx))*dx;
            end
            if (j==1)
                fs=0.5*(v_guess(cell)+v_bottom)*dx;
            end
            
            dn = Nu/dy*dx;
            ds = Nu/dy*dx;
            
            if (i<Nx)
                de=Nu/dx*dy;
            end
            if (i==Nx)
                de=2*Nu/dx*dy; % only for the top boundary cells
            end
            
            if (i>1)
                dw=Nu/dx*dy;
            end
            if (i==1)
                dw=2*Nu/dx*dy; % only for the bottom boundary cells
            end
            
            aw(cell)=dw+max(fw,0);
            ae(cell)=de+max(0,-fe);
            an(cell)=dn+max(0,-fn);
            as(cell)=ds+max(fs,0);
        
            cell_p=cell;
            
            su(cell)=(p_guess(cell_p)-p_guess(cell_p+Nx))/dy*dx*dy;
            sp(cell)=0.0;
            
            if (i==1) % left boundary
                su(cell)=su(cell)+aw(cell)*v_left;
            end
            
            if (i==Nx)
                su(cell)=su(cell)+ae(cell)*v_right;
            end
            
            if (j==1)
                su(cell)=su(cell)+as(cell)*v_bottom;
            end
            
            if (j==Ny)
                su(cell)=su(cell)+an(cell)*v_top;
            end
            
            
            ap(cell)=aw(cell)+ae(cell)+as(cell)+an(cell)-sp(cell);
            dv(cell)=1.0/dy/ap(cell); % d coefficient to be used in pressure correction equation
        end
    end
    
    % Setting up matrix
    A = zeros(N_v,N_v);
    for cell=1:N_v
        A(cell,cell)=ap(cell);
        if ((cell-1)>0)
            A(cell,cell-1)=-aw(cell);
        end
        if ((cell+1)<=N_v)
            A(cell,cell+1)=-ae(cell);
        end
        if ((cell-Nx)>0)
            A(cell,cell-Nx)=-as(cell);
        end
        if ((cell+Nx)<=N_v)
            A(cell,cell+Nx)=-an(cell);
        end
    end
    
    % Setting up RHS
    b = zeros(N_v,1);
    for cell=1:N_v
        b(cell) = su(cell);
    end
    
    % Solve for v_guess
    v_guess_updated = A\b ;
    
    % Setting up coefficients for pressure correction
    ae = zeros(N,1);
    aw = zeros(N,1);
    as = zeros(N,1);
    an = zeros(N,1);
    su = zeros(N,1);
    sp = zeros(N,1);
    
    for j=1:Ny
        for i=1:Nx
            cell=(j-1)*(Nx)+i;
            
            cell_u=cell-(j - 1);
            cell_v=cell;
            
            if (i<Nx)
                ae(cell)=du(cell_u)*dy;
                su(cell)=su(cell)-dy*u_guess_updated(cell_u); % all of the u* become source terms, read in the book
            end
            if (i==Nx)
                ae(cell)=0.0;
            end
            
            if (i>1)
                aw(cell)=du(cell_u-1)*dy;
                su(cell)=su(cell)+dy*u_guess_updated(cell_u-1);
            end
            if (i==1)
                aw(cell) = 0.0;
            end
            
            if (j>1)
                as(cell)=dv(cell_v-Nx)*dx;
                su(cell)=su(cell)+dx*v_guess_updated(cell_v-Nx);
            end
            if (j==1)
                as(cell) = 0.0;
            end
            
            if (j<Ny)
                an(cell)=dv(cell_v)*dx;
                su(cell)=su(cell)-dx*v_guess_updated(cell_v) ;
            end
            if (j==Ny)
                an(cell)=0.0;
            end
            
            sp(cell)=0;
            ap(cell)=aw(cell)+ae(cell)+as(cell)+an(cell)-sp(cell);
        
        end
    end
    
    % Setting up the matrix
    A = zeros(N,N);
    for cell=1:N
        A(cell,cell)=ap(cell);
        if ((cell-1)>0)
            A(cell,cell-1)=-aw(cell);
        end
        if ((cell+1)<=N)
            A(cell,cell+1)=-ae(cell);
        end
        if ((cell-Nx)>0)
            A(cell,cell-Nx)=-as(cell);
        end
        if ((cell+Nx)<=N)
            A(cell,cell+Nx)=-an(cell);
        end
    end
    
    
    % Setting up the RHS
    b = zeros(N,1);
    for cell=1:N
        b(cell)=su(cell);
    end
    
    % Solve for p_prime
    p_prime=A\b ;
    
    % Update p_guess
    for j=1:Ny
        for i=1:Nx
            cell=(j-1)*Nx+i;
            
            p_guess(cell)=p_guess(cell)+p_prime(cell);
            
        end
    end
    
    % Obtain u_prime and correct u_guess_updated to get u
    for j=1:Ny
        for i=1:Nx-1
            cell=(j-1)*(Nx-1)+i;
            cell_p=cell+j-1;
        
            u_prime(cell)=du(cell)*(p_prime(cell_p)-p_prime(cell_p+1));
            u_guess(cell)=u_guess_updated(cell)+u_prime(cell);
            
        end
    end
    
    % Obtain v_prime and correct v_guess_updated
    for j=1:Ny-1
        for i=1:Nx
            cell=(j-1)*Nx+i;
            cell_p=cell;
            
            v_prime(cell)=dv(cell)*(p_prime(cell_p)-p_prime(cell_p+Nx));
            v_guess(cell)=v_guess_updated(cell)+v_prime(cell);
            
        end
    end
    
    residual_p=mean(abs(p_prime)); % mean value of p_prime
    residual_u=mean(abs(u_prime)); % mean value of u_prime
    residual_v=mean(abs(v_prime)); % mean value of v_prime
    residual=max([residual_p residual_u residual_v]);
    
    n=n+1;
    disp(n); 
    
    
end

% Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% X and Y locations for u
x=zeros(Nx-1,1);
x(1)=x_min+dx;
for i=2:Nx-1
    x(i)=x(i-1)+dx;
end

y=zeros(Ny,1);
y(1)=y_min+dy/2;
for i=2:Ny
    y(i)=y(i-1)+dy;
end

[Xu,Yu]=meshgrid(x,y); % x,y locations for plotting u

figure(1) % Velocity contours
u=(reshape(u_guess,Nx-1,Ny))'; % note the transpose since MATLAB reshapes in column-major fashion
contour(Xu,Yu,u,[-4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10]);
colorbar;

figure(4) % u(0,y)/U
i=Nx/2;
for j=1:Ny
    cell=(j-1)*(Nx-1)+i;
    
    ucenter(j)=u_guess(cell)/U;
end

plot(y,ucenter)
xlabel('y (m)')
ylabel('u(0,y)/U')

% X and Y locations for v
x=zeros(Nx,1);
y=zeros(Ny-1,1);
x(1)=x_min + dx/2;
for i=2:Nx
    x(i)=x(i-1)+dx;
end

y(1)=y_min+dy;
for i=2:Ny-1
    y(i)=y(i-1)+dy;
end

[Xv,Yv]=meshgrid(x,y); % x,y locations for plotting v

figure(2) % Velocity 
v=(reshape(v_guess,Nx,Ny-1))'; % note the transpose since MATLAB reshapes in column-major fashion
contour(Xv,Yv,v,10);
colorbar;

figure(5)% v(x,0)/U
j=Ny/2;
for i=1:Nx
    cell=(j-1)*Nx+i;
    
    vcenter(i)=v_guess(cell)/U;
end

plot(x,vcenter)
xlabel('x (m)')
ylabel('v(x,0)/U')

figure(3) % Pressure 
p=(reshape(p_guess,Nx,Ny))'; % note the transpose since MATLAB reshapes in column-major fashion
contour(X,Y,p);
colorbar;
























