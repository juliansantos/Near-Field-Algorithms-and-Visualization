
load dataTest1.mat
load Aut_clean.mat

% 2-Dim DFT

clc; clear all; 
A = [1,2;3,4]; % Sample matriz: Size(M,N)
ft_field = DFT2(A)
ift_field = IDFT2(A)


function [y]=DFT2(A)
    y = zeros(size(A)); % Preallocation output
   [M,N] = size(A);  

    for p = 0:M-1
        for q = 0:N-1
            for j = 0:M-1
                for k = 0:N-1
                    y(p+1,q+1) = exp(-2*pi*1i/M)^(j*p)*exp(-2*pi*1i/N)^(k*q)*A(j+1,k+1)+y(p+1,q+1);
                end
            end
        end
    end
end

% 2-Dim IDFT

function [y]=IDFT2(A)
    y = zeros(size(A)); % Preallocation output
    [M,N] = size(A); 
    for p = 1:M
        for q = 1:N
            for j = 1:M
                for k = 1:N
                    y(p,q) = (1/(M*N))*exp(2*pi*1i/M)^((j-1)*(p-1))*exp(2*pi*1i/N)^((k-1)*(q-1))*A(j,k)+y(p,q);
                end
            end
        end
    end
end
 
