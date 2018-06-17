% -----------------------------------------------------------------  %
% Matlab Programs included the Appendix B in the book:               %
%  Xin-She Yang, Engineering Optimization: An Introduction           %
%                with Metaheuristic Applications                     %
%  Published by John Wiley & Sons, USA, July 2010                    %
%  ISBN: 978-0-470-58246-6,   Hardcover, 347 pages                   %
% -----------------------------------------------------------------  %
% Citation detail:                                                   %
% X.-S. Yang, Engineering Optimization: An Introduction with         %
% Metaheuristic Application, Wiley, USA, (2010).                     %
%                                                                    % 
% http://www.wiley.com/WileyCDA/WileyTitle/productCd-0470582464.html % 
% http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0470582464.html  %
% -----------------------------------------------------------------  %
% ===== ftp://  ===== ftp://   ===== ftp:// =======================  %
% Matlab files ftp site at Wiley                                     %
% ftp://ftp.wiley.com/public/sci_tech_med/engineering_optimization   %
% ----------------------------------------------------------------   %


% Firefly Algorithm by X-S Yang (Cambridge University  @ 2008)       %
% Usage: B5_firefly([number_of_fireflies,MaxGeneration])             %
%  eg:   B5_firefly([25,20]);                                        %
% For constrained optimization, please see: fa_constrained_demo.m    %
% -------------------------------------------------------------------%
function [best]=B5_firefly(instr)
% n=number of fireflies
% MaxGeneration=number of pseudo time steps
filename = 'elektroforize';
x = imread([filename '.jpg']);
%img=rgb2gray(x(1:200,1:200,1:3));
x=imresize(x,0.5,'bicubic');
img=rgb2gray(x(1:50,1:50,1:3));
figure(1);
image(x);
if nargin<1,   instr=[10 50];     end
n=instr(1);  MaxGeneration=instr(2);
rand('state',0);  % Reset the random generator
% ------ Four peak functions ---------------------
funstr='x+y';
% Converting to an inline function
f=vectorize(inline(funstr));
% range=[xmin xmax ymin ymax];
range=[0 200 0 200];
range=[0 50 0 50];
% ------------------------------------------------
alpha=0.5;      % Randomness 0--1 (highly random)
gamma=0.001;      % Absorption coefficient
% ------------------------------------------------
% Grid values are used for display only
Ndiv=100;
dx=(range(2)-range(1))/Ndiv; dy=(range(4)-range(3))/Ndiv;
[x,y]=meshgrid(range(1):dx:range(2),...
               range(3):dy:range(4));
z=f(x,y);
% Display the shape of the objective function
%figure(1); 
%image(img(range(1)+1:range(2),range(3)+1:range(4)));
%surfc(x,y,z);
% ------------------------------------------------
% generating the initial locations of n fireflies
[xn,yn,Lightn]=init_ffa(n,range);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized
 figure(2);
% Iterations or pseudo time marching
[Gx,Gy]=imgradientxy(img);
Gx=abs(Gx);
Gy=abs(Gy);
hist=zeros(50)
for i=1:MaxGeneration     %%%%% start iterations
    % Show the contours of the function
    %contour(x,y,z,15); hold on;
    % Evaluate new solutions
        for j=1:n
            %tikslo funkcija
            zn(j)=(Gx(xn(j),yn(j))+Gy(xn(j),yn(j)))/510;;
        end
    % Ranking the fireflies by their light intensity

    [Lightn,Index]=sort(zn);
    xn=xn(Index); yn=yn(Index);
    xo=xn;   yo=yn;    Lighto=Lightn;
    % Trace the paths of all roaming  fireflies
    
    plot(xn,yn,'.','markersize',10,'markerfacecolor','g');
    axis ij
    %pause
    % Move all fireflies to the better locations
    [xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,...
            Lighto,alpha,gamma,range);
    drawnow;
    % Use "hold on" to show the paths of fireflies
    hold off;
    for k=1:n
       hist(xn(k),yn(k))=hist(xn(k),yn(k))+1; 
    end
end   %%%%% end of iterations
figure(3)
    hold on
    %Kas per?
    imagesc(img);
    plot(xn,yn,'.','markersize',10,'markerfacecolor','b');
    axis ij
   hold off
best(:,1)=xo'; best(:,2)=yo'; best(:,3)=Lighto';
figure(4);
surf(hist);
contourf(hist,10);
% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [xn,yn,Lightn]=init_ffa(n,range)
xn=randi([range(1)+1,range(2)],1,n);
yn=randi([range(3)+1,range(4)],1,n);
Lightn=zeros(size(yn));
% Move all fireflies toward brighter ones

function [xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,...
    Lighto,alpha,gamma,range)
ni=size(yn,2); nj=size(yo,2);
for i=1:ni,
% The attractiveness parameter beta=exp(-gamma*r)
    for j=1:nj,
    %jonvabalio sviesumas
    r=sqrt((xn(i)-xo(j))^2+(yn(i)-yo(j))^2);
        if Lightn(i)<Lighto(j), % Brighter and more attractive
        beta0=0.1;     
        beta=beta0*exp(-gamma*r.^2);
        xn(i)=round(xn(i).*(1-beta)+xo(j).*beta+alpha.*(rand-0.5));
        yn(i)=round(yn(i).*(1-beta)+yo(j).*beta+alpha.*(rand-0.5));
        end
    end % end for j
end % end for i
[xn,yn]=findrange(xn,yn,range);
% Make sure the fireflies are within the range

function [xn,yn]=findrange(xn,yn,range)
for i=1:length(yn),
   if xn(i)<=range(1), xn(i)=range(1); end
   if xn(i)>=range(2), xn(i)=range(2); end
   if yn(i)<=range(3), yn(i)=range(3); end
   if yn(i)>=range(4), yn(i)=range(4); end
end
