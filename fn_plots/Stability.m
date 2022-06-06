% Compute eigenvalues of the spatial discretisation operator. The
% eigevalues, scaled by the time step dt, should lie in the stability
% region of the time-discretisation operator (RK4 in my case). In addition,
% the pseudospectra must lie within this region too. 
% 
% Problem:
% v_t + v_x = 0 => Chebyshev spatial discret. => v_t = - 2 a Dc v, with a
% the advection velocity, Dc the Chebyshev disc. matrix and v the advection
% velocity modes discretised in X = (0, 1] using N_c + 1 Chebyshev points.
% ----------------------------------------------------------------------- %
% close all

N_c     =   200;
figure
a       =   1 / 0.2;
dt      =   1E-4;

[Dc, X] =   Cheb(N_c);  
Dc      =   2 * a * Dc(2:N_c,2:N_c) * dt;    % Removing the BC
X       =   (X(2:N_c) + 1)/2;   % X \in (0,1]

[V, Lam] =   eig(Dc);
[~, ii]     =   sort(diag(Lam));
V           =   V(:,ii);
e           =   diag(Lam(ii,ii)); 

% Runge-Kutta stability region
a = 4;
x_min = -a; x_max = a; x = -a:a/100:a; 
y_min = -a; y_max = a; y  = -a:a/100:a; 
[x,y] = meshgrid(x,y);
z = x+1i*y;
R1 = 1+z;
R2 = 1+z+(1/2)*z.^2;
R3 = 1+z+(1/2)*z.^2+(1/6)*z.^3;
R4 = 1+z+(1/2)*z.^2+(1/6)*z.^3+(1/24)*z.^4;
R5 = 1+z+(1/2)*z.^2+(1/6)*z.^3+(1/24)*z.^4+(1/120)*z.^5;


figure(1); hold on; box on; axis square
title('Stability region'); xlim([x_min, x_max])
plot([x_min x_max],[0 0],'k','LineWidth',.5, 'HandleVisibility','off')
plot([0 0],[y_min y_max],'k','LineWidth',.5, 'HandleVisibility','off')
contour(x,y,abs(R1),[1 1],'DisplayName','Euler','LineWidth',1,'LineColor','r');
contour(x,y,abs(R4),[1 1],'DisplayName','RK4','LineWidth',1,'LineColor','b');
contour(x,y,abs(R5),[1 1],'DisplayName','RK5','LineWidth',1,'LineColor','g');
legend('Location','northeastoutside')
xlabel('Re($\lambda$)'); ylabel('Im($\lambda$)')
hold on;
plot3(real(e), imag(e), ones(size(e)), 'kx', 'markersize', 10)


%% Spectral radius as fnc of Nc
figure

N   =   100;
dt  =   1./1000;
a   =   - 1/0.2;

x   =   chebpts(N, [0, 1]);
x   =   x(2:N);

Dop = chebop(@(u) (u));
D = Dop(N);
D = -D(2:N,2:N);

[V,Lam] = eig(D);
Lam     = diag(Lam);


plot(real(Lam*dt), imag(Lam*dt),'r.')
eitheta = exp(1i*linspace(0,2*pi,1000));
hold on
plot(eitheta-1,'k')
axis 'equal'
xylim = dt*Lmax + 0.1;
axis([-xylim,0.1,-xylim,xylim])


%% Plot eigenvalues
figure('Units','normalized','OuterPosition',[.2,.2,.5,.8]); tiledlayout(1,1);
nexttile; 
loglog(abs(e),'o','MarkerSize',3);hold on
ylabel('eigenvalue'); title(['$N_c=$', num2str(N_c)])
legend(['D: max $|\lambda| = $',num2str(max(abs(e))/N_c,3),'$N_c$'])

nexttile; hold on
xx = -1:.01:1; 
vv = polyval(polyfit(X, [0; V(:, N_c/2-1); 0], N_c), xx); 
vv2 = polyval(polyfit(X, [0; V2(:, N_c/4-1); 0], N_c), xx); 
plot(xx, vv); plot(xx, vv2)