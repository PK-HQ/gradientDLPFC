function [X, surrogateY]=genSurrogateY(modelStruct,X,Y)
%The random function simulates new random response values,
%equal to the mean prediction plus a random disturbance 
%with the same variance as the training data.
surrogateY=random(modelStruct,X);
end