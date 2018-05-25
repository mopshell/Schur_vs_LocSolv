function [x] = boundProjectstan(x,boptions)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
LB = boptions.LB;
UB = boptions.UB;
x(x < LB) = boptions.l;
x(x > UB) = boptions.u;

x = x(:);

end

