%find eigenvalues of stress tensor 

%
% Inputs:
%   m: 8x1 vector of model parameters:
%       vertical semi-diameter (length)
%       horizontal semi-diameters (length)
%       dip (degrees)
%       strike (degrees)
%       X1 position (length)
%       X2 position (length)
%       X3 position (length)
%       strength, pressure change (force/length^2) or volume change (m^3)
%
% xyz: 3x1 vector of observation coordinates
%  nu: Poisson's ratio
%  mu: Shear modulus (force/length^2)
%

% Outputs:
%   u : 3x1 vector of displacements [ux; uy; uz];
%   D : 9x1 vector of the components of the deformation gradient tensor
%       [dux/dx; duy/dx; duz/dx; dux/dy; duy/dy; duz/dy; dux/dz; duy/dz; duz/dz]
%   E : 6x1 vector of strain tensor components
%       [e11; e12; e13; e22; e23; e33];
%   S : 6x1 vector of stress tensor components
%       [s11; s12; s13; s22; s23; s33];

M = [2000 2000 0  90 0 0 -10000 434];
xyz= [0 0 0]';
nu = 1/3;
mu = 1;

[u,D, strain, stress]=spheroid(M, xyz, nu,mu) ;

% create cauchy stress tensor
stressMat= [ stress(1) stress(2) stress(3)
             0         stress(4) stress(5)
             0         0         stress(6)];

% find eigenvectors (V) and eigenvalues (D)
[V,D]=eig(stressMat);
D=diag(D);
[minEigVal, loc]=max(D);

%when you increase pressure M(8), u(3) goes down. 

% When you apply a positive pressure change the ground surface goes down, so it's the opposite I think. 
% Negative is compression so we want the most positive value.
% % I also found that it says "The sign convention is such that a negative pressure change produces 
% an "inflationary" pattern of deformation" in the document. 
% 
% So we want the most positive value, max().

space=0:-500:-10000;

traj=eigGrid(M, n, nu, mu); 
function traj = eigGrid(m, space, nu, mu)


a = m(1);
b = m(2);
dip = mod(m(3),180);
ct = cosd(dip);
st = sind(dip);
strike = 90 - m(4);
cs = cosd(strike);
ss = sind(strike);
X1 = m(5);
X2 = m(6);
X3 = m(7);
strength = m(8);


R1 = [st 0 ct;0 1 0; ct 0 -st];
R2 = [cs -ss 0;ss cs 0; 0  0 1];
n = 50;
[X,Y,Z] = ellipsoid(0,0,0,b,b,a,n);

coords = [X(:) Y(:) Z(:)]';
epts = bsxfun(@plus,R2*R1*coords,[X1;X2;X3]);
nhat = R2*R1*bsxfun(@rdivide,2*coords,[b^2;b^2;a^2]);
nhat = bsxfun(@rdivide,nhat,sqrt(sum(nhat.^2)));

X(:) = epts(1,:)%/1000;
Y(:) = epts(2,:)%/1000;
Z(:) = epts(3,:)%/1000;
hs = surf(X,Y,Z);
alpha 0.5

axis equal
L = a/1000+1;
%axis([X1/1000+[-L L] X2/1000+[-L L] X3/1000+[-L L]])

hold on


[x,y,z ]= meshgrid( min(min(X))-1000:500:1000, min(min(Y))-1000:500:1000, min(min(Z))-1000:500:1000) ;
loc=[x(:), y(:), z(:)]';
n=max(size(loc));
[nx, ny, nz]=size(x);


xEig=zeros(n,1);
yEig=zeros(n,1);
zEig=zeros(n,1);

for i = 1 : n
    xyz=loc(:,i);
    [u,D, strain, stress]=spheroid(m, xyz, nu,mu) ;
    % create cauchy stress tensor
    stressMat= [ stress(1) stress(2) stress(3)
                 0         stress(4) stress(5)
                 0         0         stress(6)];

    % find eigenvectors (V) and eigenvalues (D)
    [V,D]=eig(stressMat);
    D=diag(D);
    [minEigVal, minloc]=max(D);
    
    
    xEig(i)=V(minloc, 1);
    yEig(i)=V(minloc, 2);
    zEig(i)=V(minloc, 3);
   
end

xEig=reshape(xEig, nx,ny, nz);
yEig=reshape(yEig, nx,ny, nz);
zEig=reshape(zEig, nx,ny, nz);


traj=[xEig, yEig, zEig];

nStream=15;
points=randi([1,size(Z,1)],1,nStream);

startx=X %(points);
starty=Y %(points);
startz=Z %(points);

vert=stream3(x,y,z, xEig, yEig, zEig, startx, starty, startz);
streamline(vert);
view(3)


end
    


