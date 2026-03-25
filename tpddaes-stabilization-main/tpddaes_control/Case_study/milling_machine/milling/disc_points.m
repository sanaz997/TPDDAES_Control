function [points]=disc_points(phi_en,phi_ex,omega,N)
points=[];
T=60/N/omega;
p1_en=60/2/pi/omega*(phi_en-1*2*pi/N);
p1_ex=60/2/pi/omega*(phi_ex-1*2*pi/N);

p2_en=60/2/pi/omega*(phi_en-2*2*pi/N+2*pi);
p2_ex=60/2/pi/omega*(phi_ex-2*2*pi/N+2*pi);

if p1_en>0 && p1_en<T
    points=[points p1_en];
end

if p1_ex>0 && p1_ex<T
    points=[points p1_ex];
end

if p2_en>0 && p2_en<T
    points=[points p2_en];
end

if p2_ex>0 && p2_ex<T
    points=[points p2_ex];
end

if ~isempty(points)
    points=sort(points);
end
end