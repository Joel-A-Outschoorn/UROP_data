clear
clc

x = 0:0.01:10;

y = x.^2;

Int = simpsons(y,0,10,[]);

X = 0:2:10;

Y = X.^2;

X_inter = 5;
Y_inter = interp1(X,Y,X_inter)

