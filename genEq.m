function F = genEq(x,M,anzahl_zweige)

c=x; % c = vector of variables, needs to match the size of inputted solutions-vector

for i = 1:size(x) % fills vector c with x(1) ... x(N)// just as an example to test if the system works..
    c(i) = x(i);
end

for i = anzahl_zweige+1:2*anzahl_zweige % the volume flow is a square function
 c(i) = x(i)^2;
end
    
F = M*c  % creates functions that have to be solved (F(x)=0)



