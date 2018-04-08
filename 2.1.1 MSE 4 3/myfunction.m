function F = myfunction(e,arfa,d,r0,gama,beta,r,theta)
F = @(fai) [sin(fai-theta)-((d+e.*sin(theta-arfa))./(r+r0.*sqrt((((cos(fai-beta)).^2)./(cos(gama)).^2)+(sin(fai-beta)).^2)))];
x0=[0];
fai=fsolve(F,x0)
end
