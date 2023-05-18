function vector = runge_kutta_4_vector(f, p0, detaT)

k1 = detaT*f(p0); 
k1 = k1/norm(k1);
p1 = p0 + k1/2; 
k2 = detaT*f(p1); 
k2 = k1/norm(k2); 
p2 = p0 + k2/2;
k3 = detaT*f(p2); 
k3 = k1/norm(k3); 
p3 = p0 + k3;
k4 = detaT*f(p3);
k4 = k1/norm(k4); 
vector = (k1 + 2*k2 + 2*k3 + k4)/6;
vector = vector/norm(vector);

end