function [BW,laikas]=jonvabaliai(filename);
tic
%filename = 'dataset\001';
x = imread([filename '.bmp']);
%img=rgb2gray(x(1:200,1:200,1:3));
%x=imresize(x,0.5,'bicubic');
img=rgb2gray(x(1:120,1:120,1:3));
%figure(1);
%image(x);
%Jonvabaliu skaicius ir iteraciju skaicius
n=20;  
MaxGeneration=100;
rand('state',0);  % Reset the random generator
%isskiriamos ribos
[nrow, ncol] = size(img);
% ------------------------------------------------
alpha=0.5;      % Randomness 0--1 (highly random)
gamma=0.001;      % Absorption coefficient
% Jonvabaliu inicilizacija
[xn,Lightn]=init_ffa(n);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized
hist=imhist(img); % Compute the histogram
N=sum(hist); % sum the values of all the histogram values
max=0; %initialize maximum to zero
%%================================================================================================
for i=1:255
    P(i)=hist(i)/N; %Computing the probability of each intensity level
end

for i=1:MaxGeneration     %%%%% Pagrindinis ciklas
    % Tikriname sprendimus
        for j=1:n
            %tikslo funkcija 
            w0=sum(P(1:xn(j))); % Probability of class 1 (separated by threshold)
            w1=sum(P(xn(j)+1:255)); %probability of class2 (separated by threshold)
            u0=dot([0:xn(j)-1],P(1:xn(j)))/w0; % class mean u0
            u1=dot([xn(j):254],P(xn(j)+1:255))/w1; % class mean u1
            zn(j)=w0*w1*((u1-u0)^2); % compute sigma(zn(j) i.e variance(between class)
        end
    % Isrusiuojame pagal tikslo funkcija
    [Lightn,Index]=sort(zn);
    xn=xn(Index);
    xo=xn; Lighto=Lightn;
    % Sekame jonvabalius
    %pause
    % Judiname jonvabalius pagal sviesos koeficientus
    [xn]=ffa_move(xn,Lightn,xo,Lighto,alpha,gamma);
end   %%%%% end of iterations
%%
%atvaizduojame vienas prie kito rezultata
[Lightn,Index]=sort(zn);
threshold=xn(1);
BW=im2bw(x,threshold/255);
%figure(1);
%imshowpair(img,BW,'montage');
% [Gx,Gy]=imgradientxy(BW);
% Gx=abs(Gx);
% Gy=abs(Gy);
% G=Gx+Gy;
% BW=int8(G>1).*255;
temp=ones(120);
BW=uint8((temp-BW).*255);
laikas=toc;
end
% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [xn,Lightn]=init_ffa(n)
xn=randi([1,255],1,n);
Lightn=zeros(size(xn));
% Move all fireflies toward brighter ones
end
%[xn]=ffa_move(xn,Lightn,xo,Lighto,alpha,gamma,range);
function [xn]=ffa_move(xn,Lightn,xo,Lighto,alpha,gamma)
ni=size(xn,2); nj=size(xo,2);
for i=1:ni,
% The attractiveness parameter beta=exp(-gamma*r)
    for j=1:nj,
    %jonvabalio sviesumas
    r=sqrt(xn(i)-xo(j))^2;
        if Lightn(i)<Lighto(j), % Brighter and more attractive
        beta0=0.1;     
        beta=beta0*exp(-gamma*r.^2);
        xn(i)=round(xn(i).*(1-beta)+xo(j).*beta+alpha.*randi([-5 5]));
        if(xn(i)>254)
            xn(i)=254;
        end
        if(xn(i)<=0)
            xn(i)=0;
        end
        end
    end % end for j
end % end for i

end
