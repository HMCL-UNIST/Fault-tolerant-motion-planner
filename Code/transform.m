
function X_out = transform(X_in)

X_out = zeros(3,1);

U =     X_in(4);
V =     X_in(5);
W =     X_in(6);

phi =   X_in(10);
theta = X_in(11);
psi =   X_in(12);

C1=[1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
C2=[cos(theta) 0 -sin(theta) ; 0 1 0; sin(theta) 0 cos(theta)];
C3=[cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
C=(C1*C2*C3)';

X_out = C*[U; V; W];
X_out(3) = -X_out(3);

