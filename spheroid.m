function varargout = spheroid(m,xyz,nu,mu,flag)
% spheroid  [u,D,E,S] = spheroid(m,xyz,nu,mu,[strength_type])
%
% Internal and surface deformation from a dipping spheroid in a half space.
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
% Optional Input:
%   strength_type: Flag to denote the strength term as either pressure or
%                  volume change, taking one of two values: 'pressure' or
%                  'volume'. If omitted, strength term defaults to
%                  pressure change.
%
% If a single input, the model vector, is given, the spheroid is plotted
% and the pressure on its internal surface is calculated.  The color scale
% depicts the departure from uniformity of the surface pressure.
%
% If no inputs are given, 15 cases with different spheroid geometries and
% observation points are tested and the results summarized.
%
% Outputs:
%   u : 3x1 vector of displacements [ux; uy; uz];
%   D : 9x1 vector of the components of the deformation gradient tensor
%       [dux/dx; duy/dx; duz/dx; dux/dy; duy/dy; duz/dy; dux/dz; duy/dz; duz/dz]
%   E : 6x1 vector of strain tensor components
%       [e11; e12; e13; e22; e23; e33];
%   S : 6x1 vector of stress tensor components
%       [s11; s12; s13; s22; s23; s33];
% Notes:
%
%  (1) [X1,X2,X3] are the coordinates of the center of the spheroid.
%      X3, therefore, should always be negative.
%
%  (2) Strike is reckoned clockwise from north (0°). The spheroid
%      dips in the strike direction, with 0° being horizontal and 90°
%      being vertical.
%
%  (3) If a = b, the solution reduces to the spherical point source (aka
%      the "Mogi" or "Anderson" source.
%
%  (4) The sign convention is such that a negative pressure change produces
%      an "inflationary" pattern of deformation.
%
%  (5) Called with a single argument, the model vector, the code plots the
%      spheroid, evaluates the pressure change on the internal cavity wall,
%      and displays the departure from uniformity with a color scale.
%
% This code corrects and extends the expressions give in "Deformation From
% Inflation of a Dipping Finite Prolate Spheroid in an Elastic Half-Space as
% a Model for Volcanic Stressing", by Xue-Min Yang. Paul Davis, and James
% Dieterich, Journal of Geophysical Research, Volume 93, Number B5, pages
% 4249 - 4257, May 10, 1998. I have written the code to correct errors
% found in the original paper, to add expressions for the deformation gradient
% tensor, from which strains and stresses are derived, and to express the
% strength term as either a pressure or volume change. I have also
% eliminated singularities where possible, by evaluating the limits as dip
% goes to 90 degrees, as x2 - X2 goes to zero, and as R1 + r3bar goes to
% zero.

% Version 1.11, October 10, 2018 USGS (PFC)

% As of 2018/10/10, the function returns only the real part of the
% calculated displacements and derivatives in cases where rounding errors
% lead to a non-zero imaginary part.

ZERO_THRESHOLD = 1e-11;

if nargin == 0
    test
    return
end

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

if nargin == 1
    varargout{1} = draw;
    return
end

c = sqrt(a^2 - b^2);
b2 = b^2;
ct2 = ct^2;
ct2 = ct^2;
C0 = 4*(1 - nu);
C1 = C0*(1 - 2*nu);
C2 = C0 - 1;
C3 = C0 + 1;
C4 = C0 - 3;

xyz(1,:) = xyz(1,:) - X1;
xyz(2,:) = xyz(2,:) - X2;
x1 = cs*xyz(1,:) + ss*xyz(2,:);
x2 = -ss*xyz(1,:) + cs*xyz(2,:);
x3 = xyz(3,:);
x3hat = x3 - X3;
x3bar = x3 + X3;

if a == b
    C5 = -a^3*strength/(4*mu);
    if nargin > 4
        if strcmpi(flag,'volume')
            C5 = -strength/(4*pi);
        end
    end

    s1 = sqrt(x1.^2 +x2.^2 + x3hat.^2);
    s2 = sqrt(x1.^2 +x2.^2 + x3bar.^2);
    s13 = s1.^3;
    s15 = s1.^5;
    s22 = s2.^2;
    s23 = s2.^3;
    s25 = s2.^5;

    u = C5*[x1.*(1./s13 + C2./s23 + 6*x3.*x3bar./s25)
            x2.*(1./s13 + C2./s23 + 6*x3.*x3bar./s25)
            x3hat./s13 + (-C4.*x3bar - 2*X3)./s23 + 6*x3.*x3bar./s25];

    if nargout> 1
        D = C5*[1./s13 - 3*x1.^2./s15 + C2./s23.*(1 - 3*x1.^2./s22) - 6*x3.*x3bar./s25.*(1 - 5*x1.^2./s22)
                3*x1.*x2.*(-1./s15 - 1./s25.*(C2 - 10*x3.*x3bar./s22))
                -3*x1.*((x3 - X3)./s15 + (2*x3 + x3bar.*(C3 - 10*x3.*x3bar./s22))./s25)
                3*x1.*x2.*(-1./s15 - 1./s25.*(C2 - 10*x3.*x3bar./s22))
                1./s13 - 3*x2.^2./s15 + C2./s23.*(1 - 3*x2.^2./s22) - 6*x3.*x3bar./s25.*(1 - 5*x2.^2./s22)
                -3*x2.*(x3hat./s15 + (2*x3 + x3bar.*(C3 - 10*x3.*x3bar./s22))./s25)
                3*x1.*(-x3hat./s15 + (2*X3 + x3bar.*(C4 + 10*x3.*x3bar./s22))./s25)
                3*x2.*(-x3hat./s15 + (2*X3 + x3bar.*(C4 + 10*x3.*x3bar./s22))./s25)
                1./s13 - 3*x3hat.^2./s15 - C4./s23 - 3.*x3bar./s25.*(6*x3 - x3bar.*(C4 + 10*x3.*x3bar./s22))];
    end
else

    L0 = log((a - c)/(a + c));
    d = 1/(b2*c*(2*b2 - C0*a^2)*L0 - a*c^2*(4*b2 + C0*a^2) + a*b^4*(1 + nu)*L0^2);
    a1 = d*b2*(2*c*(a^2 + 2*b2) + 3*a*b2*L0);
    b1 = d*(a^2*c*C0 + b2*(c + C3*(c + a*L0)));

    C5 = strength/(4*mu);
    if nargin > 4
        if strcmpi(flag,'volume')
            C5 = strength*3/(8*b2*pi*(-2*b1*c^3 - a1*(a*L0*(1 - 2*nu) + c*C3)));
        end
    end

    r2 = x1*st + x3hat*ct;
    r3 = x1*ct - x3hat*st;
    q2 = x1*st - x3bar*ct;
    q3 = -x1*ct - x3bar*st;

    if nargout > 1
        [u1,D1] = evaluate(c);
        [u2,D2] = evaluate(-c);
        D = real(C5*(D1 - D2));
    else
        u1 = evaluate(c);
        u2 = evaluate(-c);
    end
    u = real(C5*(u1 - u2));

end

varargout{1} = [cs*u(1,:) - ss*u(2,:)
                ss*u(1,:) + cs*u(2,:)
                u(3,:)];
if nargout > 1
    D = [cs^2*D(1,:) - cs*(D(2,:) + D(4,:))*ss + D(5,:)*ss^2
         cs*(cs*D(2,:) + D(1,:)*ss) - ss*(cs*D(5,:) + D(4,:)*ss)
         cs*D(3,:) - D(6,:)*ss
         cs^2*D(4,:) + cs*(D(1,:) - D(5,:))*ss - D(2,:)*ss^2
         cs^2*D(5,:) + cs*(D(2,:) + D(4,:))*ss + D(1,:)*ss^2
         cs*D(6,:) + D(3,:)*ss
         cs*D(7,:) - D(8,:)*ss
         cs*D(8,:) + D(7,:)*ss
         D(9,:)];
    varargout{2} = D;
end

if nargout > 2
    E = zeros(6,size(D,2));
    E(1,:) = D(1,:);
    E(2,:) = 1/2*(D(2,:) + D(4,:));
    E(3,:) = 1/2*(D(3,:) + D(7,:));
    E(4,:) = D(5,:);
    E(5,:) = 1/2*(D(6,:) + D(8,:));
    E(6,:) = D(9,:);
    varargout{3} = E;
end

if nargout > 3
    S = zeros(6,size(D,2));
    theta = 2*mu*nu/(1-2*nu)*(E(1,:) + E(4,:) + E(6,:));
    S(1,:) = theta + 2*mu*E(1,:);
    S(2,:) = 2*mu*E(2,:);
    S(3,:) = 2*mu*E(3,:);
    S(4,:) = theta + 2*mu*E(4,:);
    S(5,:) = 2*mu*E(5,:);
    S(6,:) = theta + 2*mu*E(6,:);
    varargout{4} = S;
end

function [U,D] = evaluate(Z)

y1 = x1 - Z*ct;
y2 = x2;
y3 = -x3hat - Z*st;
y3bar = -x3bar + Z*st;
r3bar = r3 - Z;
q3bar = q3 + Z;
rho2 = y1.^2 + y2.^2;
R1 = sqrt(rho2 + y3.^2);
R12 = R1.^2;
R13 = R1.^3;
R2 = sqrt(rho2 + y3bar.^2);
R22 = R2.^2;
R23 = R2.^3;
R24 = R2.^4;
R25 = R2.^5;

R1r3 = R1 + r3bar;
R2q3 = R2 + q3bar;
R2y3 = R2 + y3bar;

Z2 = Z^2;

I1 = abs(R1r3) <= ZERO_THRESHOLD;
I2 = rho2 == 0;
I3 = y2 == 0;

L1 = log(R1r3);
A1star = a1./(R1.*R1r3) + b1*(L1 + (r3 + Z)./R1r3);

if st == 0
    L1(I1) = 1/2*log((x1(I1) + Z).^2);
    A1star(I1)  = -a1*x1(I1)*Z/(x1(I1).^2 - Z2).^2 - b1*((x1(I1) + Z)/(x1(I1) - Z)/2 - L1(I1));
else
    L1(I1) = log(-Z + (x3(I1) - X3)/st);
    A1star(I1) = 1/2*(a1./(-Z+(x3(I1) - X3)/st).^2 + b1*(2*L1(I1) + (x3(I1) - X3 + Z*st)./(x3(I1) - X3 - Z*st)));
end

L2 = log(R2q3);
L3 = log(R2y3);
A1 = Z./R1 + L1;
A1bar = Z./R2 - L2;
A2 = R1 - r3.*L1;
A2bar = R2 - q3.*L2;
A3 = Z*r3bar./R1 + R1;
A3bar = Z*q3bar./R2 - R2;
A1starbar =  -a1./(R2.*R2q3) - b1*(L2 + (q3 - Z)./R2q3);

if ct == 0
    f1 = C1*y1.*(Z./R2y3 + L2/2 - (R2.*(3*Z + x3bar)/2 + x3bar*Z)./rho2);
    f2 = C1*y2.*(Z./R2y3 + L2/2 - (R2.*(3*Z + x3bar)/2 + x3bar*Z)./rho2);
    f1(I1 | I2) = 0;
    f2(I1 | I2) = 0;
    f3 = C1*(R2 + x3bar.*L2);
else
    phi = atan(-(ct*q2 + (1 + st)*R2q3)./(ct*y2));
    f1 = C1*(Z*y1./R2y3 - 6*st*y2.*phi/ct2 + 3/ct2*(q2.*(st*L2 - L3) + ct*(R2 - y3bar)) - 2*ct*A2bar - 2/ct*(x3bar.*L3 + q3.*L2));
    f2 = C1*(Z*y2./R2y3 - 6*q2.*phi/ct2 - 3*y2/ct2.*(L2 - st*L3) + 2*y2.*L2 - 4*x3bar/ct.*phi);
    f3 = C1*(-2*y2.*phi/ct + 2*st*A2bar + q3.*L3 - Z + q2/ct.*(L2 - st*L3));
end

F1 = -2*x3.*(C3*st - (X3 - Z*(C3*st - (X3 - st*Z)*(q3bar./R2 + 1)./R2))./R2)./R2q3 - C2*A1bar;
F2 = -2*x3.*(A1bar*C3*st - (X3 + q3bar*Z*(X3 - st*Z)./R22)./R2) - C2*A3bar;
F1star = -2*x3.*(st*(2*b1*(R2 + Z)./R2q3 - a1./R22) + ct*q2.*(b1*(R2 + 2*Z) - a1*(q3bar./R2 + 2)./R2)./R2q3.^2)./R2;
F2star = 2*x3.*(a1*y3bar./R23 - 2*b1*(A1bar*st + ct*q2.*(R2 + Z)./R2./R2q3));
Bstar = 2*A2*b1 + a1./R1 + C2*(2*A2bar*b1 + a1./R2);
B1 = -2*C2*(X3*(Z./R2 - L2) + st*(A2bar - Z2./R2));
B2 = A3 + C0*(A2 + A2bar);
B3 = A1star.*r2 + C2*q2.*A1starbar;
B4 = 2*ct*x3.*A1starbar + F1star.*q2;

U = [b2*(st*(B3 - B4) + ct*(Bstar + F2star)) + a1*(st*(q2.*F1 - A1.*r2) - ct*(F2 + B2) + f1)
     b2*y2.*(A1star + C2*A1starbar - F1star) - a1*(y2.*(A1 - F1) - f2)
     b2*(ct*(B3 + B4) - st*(Bstar - F2star)) + a1*(-ct*(q2.*F1 + A1.*r2) - st*(F2 - B2) - B1 + f3)];

if nargout > 1

    L1dx1 = (ct + y1./R1)./R1r3;
    L1dx2 = y2./(R1.*R1r3);
    L1dx3 =  -(st + y3./R1)./R1r3;
    A1stardx1 = (-(a1./R1).*(y1./R12 + L1dx1) + b1*(2*ct + y1./R1 - (r3bar + 2*Z).*L1dx1))./R1r3;
    A1stardx2 = y2.*(-a1*(1./R1 + 1./R1r3)./R12 + b1*(1 - 2*Z./R1)./R1r3)./R1r3;
    A1stardx3 = (a1*(-L1dx3 + y3./R12)./R1 - b1*(2*st + y3./R1 + L1dx3.*(r3bar + 2*Z)))./R1r3;

    if st == 0
        L1dx1(I1) = 1./sqrt((Z - x1(I1)).^2);
        L1dx2(I1) = 0;
        L1dx3(I1) = 0;
%         A1stardx1(I1) = (a1*R13(I1)./(x1(I1).^2-Z2).^2+b1*(Z*(x1(I1).^2+Z2)./(x1(I1).^2-Z2)-R1(I1)))./(x1(I1).^2-Z2);
%         A1stardx2(I1) = 0;
%         A1stardx3(I1) = 0;
    else
        L1dx1(I1) = ct./(Z + (x3(I1) - X3)/st);
        L1dx2(I1) = 0;
        L1dx3(I1) = -st./(Z + (x3(I1) - X3)/st);
%         A1stardx1(I1) = ct*(-Z + (x3(I1) - X3)/st).^(-2).*(a1./(-Z + (x3(I1) - X3)/st) - b1*(-2*Z + (x3(I1) - X3)/st));
%         A1stardx2(I1) = 0;
%         A1stardx3(I1) = -st*(-Z + (x3(I1) - X3)/st).^(-2).*(a1./(-Z + (x3(I1) - X3)/st) - b1*(-2*Z + (x3(I1) - X3)/st));
    end
    A1stardx1(I1) = 0;
    A1stardx2(I1) = 0;
    A1stardx3(I1) = 0;
    L2dx1 = (-ct + y1./R2)./R2q3;
    L2dx2 = y2./(R2.*R2q3);
    L2dx3 =  -(st + y3bar./R2)./R2q3;
    L3dx1 = y1./(R2.*R2y3);
    L3dx2 = y2./(R2.*R2y3);
    L3dx3 =  -(1 + y3bar./R2)./R2y3;
    A1starbardx1 = a1./(R2.*R2q3).*(L2dx1 + y1./R22) - b1*(-ct./R2q3 + L2dx1.*(1 - (q3bar - 2*Z)./R2q3));
    A1starbardx2 = y2./(R2.*R2q3).*(a1./R2.*(1./R2 + 1./R2q3) - b1*(1 - (q3bar - 2*Z)./R2q3));
    A1starbardx3 = a1./(R2.*R2q3).*(L2dx3 - y3bar./R22) - b1*(-st./R2q3 + L2dx3.*(1 - (q3bar - 2*Z)./R2q3));
    A1dx1 = -Z*y1./R13 + L1dx1;
    A1dx2 = -Z*y2./R13 + L1dx2;
    A1dx3 = Z*y3./R13 + L1dx3;
    A2dx1 = y1./R1 - r3.*L1dx1 - ct*L1;
    A2dx2 = y2./R1 - r3.*L1dx2;
    A2dx3 = -y3./R1 - r3.*L1dx3 + st*L1;
    A2bardx1 = y1./R2 - q3.*L2dx1 + ct*L2;
    A2bardx2 = y2./R2 - q3.*L2dx2;
    A2bardx3 = -y3bar./R2 - q3.*L2dx3 + st*L2;
    A3dx1 = (y1 + Z*(ct - y1.*r3bar./R12))./R1;
    A3dx2 = y2.*(1 - Z*r3bar./R12)./R1;
    A3dx3 = (-y3 + Z*(-st + y3.*r3bar./R12))./R1;
    Bstardx1 = -a1*y1.*(1./R13 + C2./R23) + 2*b1*(A2dx1 + C2*A2bardx1);
    Bstardx2 = -a1*y2.*(1./R13 + C2./R23) + 2*b1*(A2dx2 + C2*A2bardx2);
    Bstardx3 = a1*(y3./R13 + C2*y3bar./R23) + 2*b1*(A2dx3 + C2*A2bardx3);

    F1dx1 = C2*(L2dx1 + y1*Z./R23) + 2*x3.*(-3*y1*Z*(X3 - st*Z)./R25 - (X3./R2.*(L2dx1 + y1./R22) - ...
            C3*st*(L2dx1.*(1 + Z./R2) + y1*Z./R23))./R2q3);
    F1dx2 = x2./R2.*(C2*(1./R2q3 + Z./R22) + 2*x3.*(-3*Z*(X3 - st*Z)./R2.^4 - X3./R2.*(1./R2q3 + 1./R2)./R2q3 + ...
            C3*st./R2q3.*(1./R2q3 + Z./R2.*(1./R2q3 + 1./R2))));
    F1dx3 = C2*(L2dx3 - y3bar*Z./R23) + 2*(Z*(X3 - st*Z)./R23 + X3./R2./R2q3 - C3*st./R2q3.*(1 + Z./R2)) + ...
            2*x3.*(3*y3bar*Z*(X3 - st*Z)./R25 - X3./R2.*(L2dx3 - y3bar./R22)./R2q3 + ...
            C3*st./R2q3.*(L2dx3.*(1 + Z./R2) - y3bar*Z./R23));

    F2dx1 = C2./R2.*(y1 + ct*Z+q3bar.*y1*Z./R22) + 2*x3.*(C3*st*(L2dx1 + y1*Z./R23) - ...
            (y1*X3 + Z*(X3 - st*Z)*(3*q3bar.*y1./R22 + ct))./R23);
    F2dx2 = x2.*(C2./R2.*(1 + q3bar*Z./R22) + 2*x3.*(C3*st./R2.*(1./R2q3 + Z./R22) - ...
            (X3 + 3*q3bar*Z*(X3 - st*Z)./R22)./R23));
    F2dx3 = C2*(st*Z - y3bar.*(1 + q3bar*Z./R22))./R2 + 2*x3.*(C3*st*(L2dx3 - y3bar*Z./R23) + ...
            (y3bar*X3 + Z*(X3 - st*Z)*(3*q3bar.*y3bar./R22 - st))./R23) + ...
            2*(X3./R2 + q3bar*Z*(X3 - st*Z)./R23 - C3*st*(Z./R2 - L2));

    if ct == 0

        F1stardx1 = -2*x3.*y1./R22.*(3*a1./R23 + 2*b1./R2y3.*(1 - (Z + R2).*(1./R2 + 1./R2y3)));
        F1stardx2 = -2*x3.*y2./R22.*(3*a1./R23 + 2*b1./R2y3.*(1 - (Z + R2).*(1./R2 + 1./R2y3)));
        F1stardx3 = 2./R2.*(x3.*(3*a1*y3bar./R24 + 2*b1./R2y3.*(y3bar./R2 + (Z + R2).*(L3dx3 - y3bar./R22))) + ...
                    a1./R22 - 2*b1*(Z + R2)./R2y3);
        F2stardx1 = -2*y1.*x3./R2.*(3*a1*y3bar./R24 - 2*b1*(Z./R22 + 1./R2y3));
        F2stardx2 = -2*y2.*x3./R2.*(3*a1*y3bar./R24 - 2*b1*(Z./R22 + 1./R2y3));
        F2stardx3 = 2*(a1./R23.*(y3bar - x3.*(1 - 3*y3bar.^2./R22)) + 2*b1*(-Z./R2 + x3.*(-Z*y3bar./R23 + L3dx3) + L2));

        f1dx1 = C1*((L2 + y1.*L2dx1)/2 - Z./R2y3.*(y1.*L3dx1 - 1) + ((x3bar*Z + R2/2.*(x3bar + 3*Z)).*(2*y1.^2./rho2 - 1) - ...
                ((x3bar + 3*Z).*y1.^2)./(2*R2))./rho2);
        f1dx1(I1) = C1/2*((Z*(R2(I1) - 2*(x3(I1) + X3)) + (x3(I1) + X3).*sqrt((x3(I1) + X3 + Z).^2))./(R2q3(I1).*(x3(I1) + ...
                    X3 + Z - sqrt((x3(I1) + X3 + Z).^2))) + L2(I1));
        f1dx1(~I1 & I2) = 0;
        f1dx2 = C1*y1.*y2.*((2*x3bar*Z + R2.*(x3bar + 3*Z))./rho2.^2 - ((x3bar + 3*Z)./(2*rho2) + Z./R2y3.^2 - 1./(2*R2q3))./R2);
        f1dx2(I1) = 0;
        f1dx2(~I1 & I2) = 0;
        f1dx3 = C1*y1.*(((((x3bar + 3*Z).*y3bar)./R2 - R2)/2 - Z)./rho2 + L2dx3/2 - Z*L3dx3./R2y3);
        f1dx3(I1) = 0;
        f1dx3(~I1 & I2) = 0;

        f2dx1 = C1*y1.*y2.*((2*x3bar*Z + R2.*(x3bar + 3*Z))./rho2.^2 - ((x3bar + 3*Z)./(2*rho2) + Z./R2y3.^2 - 1./(2*R2q3))./R2);
        f2dx1(I1) = 0;
        f2dx1(~I1 & I2) = 0;
        f2dx2 = C1*((L2 + y2.*L2dx2)/2 - Z./R2y3.*(y2.*L3dx2 - 1) + ((x3bar*Z + R2/2.*(x3bar + 3*Z)).*(2*y2.^2./rho2 - 1) - ...
                ((x3bar + 3*Z).*y2.^2)./(2*R2))./rho2);
        f2dx2(I1) = f1dx1(I1);
        f2dx2(~I1 & I2) = 0;
        f2dx3 = C1*y2.*(1./rho2.*((((x3bar + 3*Z).*y3bar)./R2 - R2)/2 - Z) + L2dx3/2 - Z*L3dx3./R2y3);
        f2dx3(I1) = 0;
        f2dx3(~I1 & I2) = 0;

        f3dx1 = C1*(y1./R2 + x3bar.*L2dx1);
        f3dx2 = C1*y2./R2.*(1 + x3bar./R2q3);
        f3dx3 = C1*(L2 + x3bar.*L2dx3 - y3bar./R2);
    else
        F1stardx1 = -2*x3.*(a1*(3*st*y1./R25 + ct./(R22.*R2q3).*((1./R2 + 1./R2q3).*(q2.*(L2dx1 + 2*y1./R22) - st) + ...
                    q2.*(L2dx1./R2q3 + y1./R23))) - b1./R2q3.*(2*(st + ct*q2./R2q3).*(y1*Z./R23 + L2dx1.*(1 + Z./R2)) - ...
                    ct./R2q3.*(st - 2*(L2dx1.*q2 - st)*Z./R2)));
        F1stardx2 = -2*x3.*y2./R2.*(a1./R22.*(3*st./R22 + ct*q2./R2q3.*(3./R22 + 1./R2q3.*(3./R2 + 2./R2q3))) - ...
                    2*b1./R2q3.*((st + ct*q2./R2q3).*(1./R2q3 + Z./R2.*(1./R2 + 1./R2q3)) + ct*q2*Z./(R2.*R2q3.^2)));
        F1stardx3 = -2*(a1./R22.*(ct./R2q3.*((1./R2 + 1./R2q3).*(x3.*(ct + q2.*(L2dx3 - y3bar./R22)) - q2) - ...
                    x3.*q2.*(2*y3bar./R23 - 1./R2q3.*(L2dx3 - y3bar./R22))) - st./R2.*(1 + 3*x3.*y3bar./R22)) + ...
                    b1./R2q3.*(st + (st + ct*q2./R2q3).*(1 + 2*Z./R2 + 2*x3.*(Z*y3bar./R23 - L2dx3.*(1 + Z./R2))) - ...
                    ct*x3./R2q3.*(ct + 2*Z./R2.*(L2dx3.*q2 + ct))));
        F2stardx1 = -2*x3.*(3*a1*y3bar.*y1./R25 - 2*b1*(st*(Z*y1./R23 + L2dx1) - ct./(R2.*R2q3).*(y1.*q2./R2 + ...
                    (R2 + Z).*(st - q2.*(y1./R22 + L2dx1)))));
        F2stardx2 = -2*x3.*y2./R2.*(3*a1*y3bar./R24 + 2*b1*(ct*q2./(R2.*R2q3).*(1 - (R2 + Z).*(1./R2 + 1./R2q3)) - ...
                    st*(1./R2q3 + Z./R22)));
        F2stardx3 = -2*(-a1./R23.*(y3bar - x3.*(1 - 3*y3bar.^2./R22)) - 2*b1*(-(ct*(q2 - x3.*(ct + q2.*L2dx3)))./R2q3 + ...
                    st*(L2 - Z./R2 + x3.*(L2dx3 - Z*y3bar./R23)) - Z*(ct*q2./(R2.*R2q3).*(1 - x3.*(L2dx3 - y3bar./R22 + ...
                    ct./q2)))));

        phidx1 = -(ct*st + (1 + st)*L2dx1.*R2q3)./(y2*ct + (ct*q2 + (1 + st)*R2q3).^2./(ct*y2));
        phidx2 = (ct*q2.*R2 + (1 + st)*(R2.*R2q3 - y2.^2))./(R2.*((ct*q2 + (1 + st)*R2q3).^2 + ct2*y2.^2)/ct);
        phidx3 = (ct2 - (1 + st)*L2dx3.*R2q3)./(y2*ct + (ct*q2 + (1 + st)*R2q3).^2./(ct*y2));
        phidx1(I3) = 0;
        phidx2(I3) = -ct./(ct*x1(I3)  - (1 + st)*(R2(I3)  - x3(I3)  - X3 + Z));
        phidx3(I3) = 0;

        f1dx1 = C1*(Z./R2y3.*(1 - y1.^2./(R2.*R2y3)) - 2*ct*(y1./R2 + ct*L2 - q3.*L2dx1) + (3*q2.*(st*L2dx1 - y1./(R2.*R2y3)) + ...
                3*st*(L2*st - L3) + ct*(3*y1./R2 - 2*(x3bar.*y1./(R2.*R2y3) - ct*L2 + q3.*L2dx1)) - 6*st*y2.*phidx1)/ct2);
        f1dx2 = C1*(-y2./R2.*(y1*Z./R2y3.^2 + 2*ct*(1 - q3./R2q3)) + (y2./R2.*(3*q2.*(st./R2q3 - 1./R2y3) + ...
                ct*(3 - 2*(x3bar./R2y3 + q3./R2q3))) - 6*st*(phi + y2.*phidx2))/ct2);
        f1dx3 = C1*(2*ct*(y3bar./R2 - L2*st + L2dx3.*q3) - (y1.*L3dx3*Z./R2y3) + (3*ct*(L3 - L2*st) - 3*q2.*(L3dx3 - st*L2dx3) + ...
                ct*(3*(1 - y3bar./R2) - 2*(x3bar.*L3dx3 + L3 - L2*st + q3.*L2dx3)) - 6*y2*st.*phidx3)/ct2);
        f2dx1 = -C1*(y2.*(y1*Z./(R2.*R2y3.^2) - 2*L2dx1) + (3*y2.*(L2dx1 - y1*st./(R2.*R2y3)) + 6*st*phi + ...
                2*(2*ct*x3bar + 3*q2).*phidx1)/ct2);
        f2dx2 = -C1*(-y2.^2./R2.*(2./R2q3-Z./R2y3.^2)-Z./R2y3-2*L2+(3*y2.^2./R2.*(1./R2q3-st./R2y3)+3*(L2-L3*st)+ ...
                2*(2*ct*x3bar+3*q2).*phidx2)/ct2);
        f2dx3 = -C1*(y2.*L3dx3*Z./R2y3 - 2*y2.*L2dx3 + (3*y2.*(L2dx3 - st*L3dx3) - 2*ct*phi + 2*(2*ct*x3bar + 3*q2).*phidx3)/ct2);
        f3dx1 = C1*(-ct*L3 + y1.*q3./(R2.*R2y3) + 2*st*(y1./R2 + ct*L2 - L2dx1.*q3) + (q2.*(L2dx1 - y1*st./(R2.*R2y3)) + ...
                st*(L2 - L3*st) - 2*y2.*phidx1)/ct);
        f3dx2 = C1*(y2./R2.*(q3./R2y3 + 2*st*(1 - q3./R2q3)) + (q2.*y2./R2.*(1./R2q3 - st./R2y3) - 2*(phi + y2.*phidx2))/ct);
        f3dx3 = C1*(q3.*L3dx3 - L3*st + 2*st*(-y3bar./R2 + L2*st - L2dx3.*q3) + (q2.*(L2dx3 - L3dx3*st) - ct*(L2 - L3*st) - ...
                2*y2.*phidx3)/ct);
    end

    B1dx1 = 2*C2*(Z*y1*(X3 - Z*st)./R23 + X3*L2dx1 - st*(ct*L2 + y1./R2 - q3.*L2dx1));
    B1dx2 = 2*C2*y2./R2.*(Z*(X3 - Z*st)./R22 + X3./R2q3 - st*(1 - q3./R2q3));
    B1dx3 = 2*C2*(-Z*(X3 - Z*st)*y3bar./R23 + X3*L2dx3 - st*(st*L2 - y3bar./R2 - q3.*L2dx3));
    B2dx1 = C0*(A2dx1 + A2bardx1) + A3dx1;
    B2dx2 = C0*(A2dx2 + A2bardx2) + A3dx2;
    B2dx3 = C0*(A2dx3 + A2bardx3) + A3dx3;
    B3dx1 = A1stardx1.*r2 + A1star*st + C2*(A1starbardx1.*q2 + A1starbar*st);
    B3dx2 = C2*A1starbardx2.*q2 + A1stardx2.*r2;
    B3dx3 = A1stardx3.*r2 + A1star*ct + C2*(A1starbardx3.*q2 - A1starbar*ct);
    B4dx1 = 2*ct*x3.*A1starbardx1 + q2.*F1stardx1 + st*F1star;
    B4dx2 = 2*ct*x3.*A1starbardx2 + q2.*F1stardx2;
    B4dx3 = 2*ct*(A1starbar + x3.*A1starbardx3) + q2.*F1stardx3 - ct*F1star;

    U1dx1 = b2*(ct*(Bstardx1 + F2stardx1) + (B3dx1 - B4dx1)*st) + a1*(f1dx1 - ct*(B2dx1 + F2dx1) + ...
            st*(F1dx1.*q2 - A1dx1.*r2 - st*(A1 - F1)));
    U1dx2 = b2*(ct*(Bstardx2 + F2stardx2) + (B3dx2 - B4dx2)*st) + a1*(f1dx2 - ct*(B2dx2 + F2dx2) + ...
            (F1dx2.*q2 - A1dx2.*r2)*st);
    U1dx3 = b2*(ct*(Bstardx3 + F2stardx3) + (B3dx3 - B4dx3)*st) + a1*(f1dx3 - ct*(B2dx3 + F2dx3) + ...
            st*(F1dx3.*q2 - A1dx3.*r2 - ct*(A1 + F1)));

    U2dx1 = b2*(A1stardx1 + C2*A1starbardx1 - F1stardx1).*y2 + a1*(f2dx1 - (A1dx1 - F1dx1).*y2);
    U2dx2 = b2*(A1star + A1starbar*C2 - F1star + (A1stardx2 + A1starbardx2*C2 - F1stardx2).*y2) + ...
            a1*(-A1 + F1 + f2dx2 - y2.*(A1dx2-F1dx2));
    U2dx3 = b2*(A1stardx3 + C2*A1starbardx3 - F1stardx3).*y2 + a1*(f2dx3 - (A1dx3 - F1dx3).*y2);

    U3dx1 = b2*((B3dx1 + B4dx1)*ct - (Bstardx1 - F2stardx1)*st) - a1*(B1dx1 - f3dx1 + ...
            ct*(st*(A1 + F1) + F1dx1.*q2 + A1dx1.*r2) - st*(B2dx1 - F2dx1));
    U3dx2 = b2*((B3dx2 + B4dx2)*ct - (Bstardx2 - F2stardx2)*st) - a1*(B1dx2 - f3dx2 + ...
            ct*(F1dx2.*q2 + A1dx2.*r2) - st*(B2dx2-F2dx2));
    U3dx3 = b2*((B3dx3 + B4dx3)*ct - (Bstardx3 - F2stardx3)*st) + a1*(f3dx3 - B1dx3 - ...
            ct*(ct*(A1 - F1) + F1dx3.*q2 + A1dx3.*r2) + st*(B2dx3 - F2dx3));

    D = [U1dx1;U2dx1;U3dx1;U1dx2;U2dx2;U3dx2;U1dx3;U2dx3;U3dx3];

end
end

function sigma = draw

    FONTNAME = 'Consolas';
    FONTSIZE = 12;
    R1 = [st 0 ct;0 1 0; ct 0 -st];
    R2 = [cs -ss 0;ss cs 0; 0  0 1];

    n = 50;
    [X,Y,Z] = ellipsoid(0,0,0,b,b,a,n);

    coords = [X(:) Y(:) Z(:)]';
    epts = bsxfun(@plus,R2*R1*coords,[X1;X2;X3]);
    nhat = R2*R1*bsxfun(@rdivide,2*coords,[b^2;b^2;a^2]);
    nhat = bsxfun(@rdivide,nhat,sqrt(sum(nhat.^2)));

    X(:) = epts(1,:)/1000;
    Y(:) = epts(2,:)/1000;
    Z(:) = epts(3,:)/1000;
    hs = surf(X,Y,Z);
    [~,G] = spheroid(m,epts,1/4,1);
    sigma = 2*nhat(1,:).*(nhat(2,:).*(G(4,:) + G(2,:)) + nhat(3,:).*(G(7,:) + G(3,:))) + ...
            2*nhat(2,:).*nhat(3,:).*(G(8,:) + G(6,:)) + nhat(1,:).^2.*(3*G(1,:) + G(5,:) + ...
            G(9,:)) + nhat(2,:).^2.*(G(1,:) + 3*G(5,:) + G(9,:)) + nhat(3,:).^2.*(G(1,:) + G(5,:) + 3*G(9,:));

    CData = zeros(size(X));
    CData(:) = (strength - sigma)/strength*100;

    set(hs,'edgecolor','none','CData',CData,'facecolor','interp')
    axis equal
    set(gca,'FontSize',FONTSIZE,'FontName',FONTNAME)
    L = a/1000+1;
    axis([X1/1000+[-L L] X2/1000+[-L L] X3/1000+[-L L]])
    xlabel('Easting (km)','FontSize',FONTSIZE,'FontName',FONTNAME)
    ylabel('Northing (km)','FontSize',FONTSIZE,'FontName',FONTNAME)
    zlabel('Vertical (km)','FontSize',FONTSIZE,'FontName',FONTNAME)
    title(sprintf('dip = %0.1f°, strike = %0.1f°, Maximum deviation: %.2f%%, Mean deviation: %.2f%%\n', ...
                   dip,strike,max(abs(CData(:))),mean(CData(:))),'FontSize',FONTSIZE,'FontName',FONTNAME);

    acb = colorbar;
    caxis([0 2])
    ylabel(acb,'Deviation from Uniform Pressue (%)','FontSize',FONTSIZE,'FontName',FONTNAME)
    set(acb,'FontSize',FONTSIZE,'FontName',FONTNAME)
    view(3)
    grid on

end
end

function test

u0 = [
    4965.5242315718669   24008.183123083669   7459.8656139218037 -20882.907513164533 -20854.669741963928  -32443.745343960363 -14173.0852079442560 -44243.7822238061560 171070.918242265560 386981.107133299290  33.657308049112935      0.000000000000       0.000000000000      0.000000000000      0.000000000000   4422.4189054516446
  -39352.9315082081810  -50524.354095122930      0.0000000000000      0.000000000000  -1311.614449180121       0.000000000000  -1010.5917910196965  -1774.0522964161712      0.000000000000      0.000000000000 -79.794375253182380      0.000000000000       0.000000000000      0.000000000000      0.000000000000  30888.3012009285080
  -43539.0587956483070 -104177.300567389790 -17654.0420245818640 -47765.263062045087 -47723.266843127734 -127641.086310222350 -28564.4794970734370 -46310.4959169799910 -12546.312596100213 -16042.436789624553 -95.970278082477307 -22599.412825045478  -80174.178160942203 -40267.722640112508 -41665.789244705156 -39761.4898522974250
     ];
D0 = [
     -10.6254632843368630 -13.6718102488900670 -0.9912987285358267 -2.1802544682981635 -2.1833301376496572  -2.9552930212368249  2.5172450325833449  5.466899640192978  249.0832385279988200 184.5058562492133700 -0.0222642612637464 -2.366477367893272 -7.649750195064291 -4.3502820952804298 -3.9276663795258973 -3.7616157510013242
      -4.9535106378531859  -4.0554535532422413  0.0000000000000000  0.0000000000000000  0.2968492959373341   0.0000000000000000  0.3413958488234479  0.455063859972705    0.0000000000000000   0.0000000000000000 -0.0072215793063832  0.000000000000000  0.000000000000000  0.0000000000000000  0.0000000000000000  0.5355274398570179
      -5.3438053327770598  -5.2995423294196220 -1.9717895807974593  7.0348903040142314  7.0227966820958523  26.9365857417870810  7.5450911060338397  9.274299426472810   -0.4792238538580499  -0.2796381903054455 -0.0030521209033807  0.000000000000000  0.000000000000000  0.0000000000000000  0.0000000000000000  0.6512380673470573
      -1.6655295694045971 -10.0777019456463480  0.0000000000000000  0.0000000000000000  0.2968492959373341   0.0000000000000000  0.1489534230055995  0.805956327888107    0.0000000000000000   0.0000000000000000 -0.0151843368235131  0.000000000000000  0.000000000000000  0.0000000000000000  0.0000000000000000  0.5355274398570179
       5.0946905642781575   5.0469368846713136 -1.8420326429749321 -6.9125811033315241 -6.8845641760131855 -10.7394059397419280 -2.9336754631818414 -5.139884793846173 -140.1358750342518800 -75.5656709141820780  0.0087416517247752 -2.366477367893272 -7.649750195064291 -3.9072128955991894 -4.2211248359826810 -0.0979083811940896
      12.1320507781120510  34.4302619230453080  0.0000000000000000  0.0000000000000000  0.4416853259179781   0.0000000000000000  0.3110669386406736  0.639875685442537    0.0000000000000000   0.0000000000000000  0.0298627634489056  0.000000000000000  0.000000000000000  0.0000000000000000  0.0000000000000000  4.5485599640796890
      -3.5284886027429345  -7.8550520363567040  1.9717895807974593  8.4345994180920592  8.4091232254470221   7.6093031293301125  0.3924021394016214  7.102543894335888   -0.0459537371364620   0.0295019306628666 -0.0152052478912535  0.000000000000000  0.000000000000000  0.0000000000000000 -0.0000000000000000 -0.8725695330793253
      19.1586585691957400  18.5614887724329730  0.0000000000000000  0.0000000000000000  0.5288756745564165   0.0000000000000000  0.0649840852176032  0.279034426124058    0.0000000000000000   0.0000000000000000  0.0421806760512103  0.000000000000000  0.000000000000000  0.0000000000000000  0.0000000000000000 -6.0944454002950739
       4.4009186440451069   9.4673528300195695  1.4166656857553794  5.9784003931414773  5.9589737843522164  10.1784958126522210 -1.0212053223420909 -1.115151169722307 -141.2438129510397600 -76.8914392518668710  0.0098596837343591  2.366477367893272  7.649750195064291  4.1287474954398089  4.0743956077542895  1.9222224262374250
     ];

nu = 1/3;
mu = 1;


M = [3000 2000 67  90 0 0 -10000 434
     2000 3000 67  90 0 0 -10000 434
     3000 2000 67  90 0 0 -10000 434
     3000 2000 90  90 0 0  -7500 434
     3000 2000 90  90 0 0  -7500 434
     2000 3000 90  90 0 0  -7500 434
     3000 2000  0  90 0 0  -7500 434
     2000 3000  0  90 0 0  -7500 434
     3000 2000  0  90 0 0 -10000 434
     2000 3000  0  90 0 0 -10000 434
     3000 2000 67 225 0 0 -10000 1
     3000 2000 90  90 0 0 -10000 434
     2000 3000 90  90 0 0 -10000 434
     2000 3000 0  90 0 0 -10000 434
     3000 2000 0  90 0 0 -10000 434
     3000 3000 0  90 0 0 -10000 434]';

xyz = [-888 3337 -6416
       -888 3337 -6416
       M(7,3)/tand(M(3,3)) 0 0
       3021 0 -2416
       3021 190 -2416
       3021 0 -2416
       6511 343 -3416
       6511 343 -3416
      -3000 0 -10000
      -2000 0 -10000
      -888 3337 -6416
         0    0     0
         0    0     0
         0    0     0
         0    0     0
     -1031 -7201 -407]';

line = repmat('-',75,1);
n = size(M,2);
check = zeros(n,1);
for i = 1 : size(M,2)
    [u,D] = spheroid(M(:,i),xyz(:,i),nu,mu);
    fprintf('Case %d, Displacements:\n',i)
    fprintf('%s\n',line)
    fprintf('%25.16f%25.16f%25.16f\n',[u,u0(:,i),u-u0(:,i)]');
    fprintf('%s\n',line)
    fprintf('Case %d, Derivatives:\n',i)
    fprintf('%s\n',line)
    fprintf('%25.16f%25.16f%25.16f\n',[D,D0(:,i),D-D0(:,i)]');
    fprintf('%s\n',line)
    check(i) = norm([u;D]-[u0(:,i);D0(:,i)]);
    if check(i) < 5e-8
        fprintf('PASS%71.16f\n',check(i));
    else
        fprintf(2,'FAIL    FAIL    FAIL    FAIL    FAIL    FAIL    FAIL     %18.16f\n',check(i));
    end
    fprintf('%s\n',line)
end
fprintf('\nSummary:\n')
fprintf('%s\n',repmat('-',34,1))
for i = 1:n
    fprintf('%3d: %18.16e',i,check(i));
    if check(i) < 5e-8
        fprintf('\t<strong>PASS</strong>\n');
    else
        fprintf(2,'\t<strong>FAIL</strong>\n');
    end
end
fprintf('%s\n',repmat('-',34,1))
end