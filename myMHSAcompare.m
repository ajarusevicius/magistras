function [ans,r,laikas]=myMHSAcompare(filename)
%filename = 'elektroforize';
% Bylos nuskaitymas
x = imread([filename '.bmp']);
test = imresize(x, 0.25);
% Vertimas á pilkus atspalvius
img = rgb2gray(test);

% Darbinio kvadrato dydis (nelyginis sk.)
kvad = 30; kvad2 = kvad^2;
%%
epclon = 0.8;
lemta1 = 0.05;
lemta2 = 0.05;
lemta3 = 0.8;
lemta4 = 0.8;
% Visø uodø svoris
c = zeros(kvad2, kvad2);
%??? Pilkos visø uodø vertës
r = zeros(kvad2, kvad2);
x = zeros(kvad2, kvad2);
% Visas vaizdas
% figure(1); imshow(img);

%%
ximg = zeros(kvad);
yimg = zeros(kvad);

% Konktreti vieta paveiksle
%y_0 = 105; x_0 = 120;

% Darbinis vaizdas
img_0 = double(img(1:kvad, 1:kvad));


%yimg1 = diff(fliplr(img_0));                           % iðvestinë y-aðimi
%ximg1 = diff(fliplr(img_0)')';                         % iðvestinë x-aðimi
%vimg1 = sqrt(ximg2(1:kvad,:).^2 + yimg2(:,1:kvad).^2); % pilna iðvestinë

%yimg2 = diff(img_0);                                   % iðvestinë y-aðimi
%ximg2 = diff(img_0')';                                 % iðvestinë x-aðimi
%vimg2 = sqrt(ximg2(1:kvad,:).^2 + yimg2(:,1:kvad).^2); % pilna iðvestinë

%% 2064-1-5467-1-10-20171227.pdf
%  [Gx,Gy] = imgradientxy(img_0);%,'CentralDifference');

% Gmag :: sqrt(Gx.^2+Gy.^2)  
% Gdir :: atan(Gy./Gx)
[Gmag, Gdir] = imgradient(img_0);                        
%figure(1);surf(Gmag); 
%figure(2);surf(img_0);

%% STUFF
%vimg2 = Gmag(1:kvad,1:kvad);
%vimg2 = vimg2./max(max(vimg2));                        % normavimas iki 1
%img_1 = imgaussfilt(img_0,2);
%img_1 = img_1-min(min(img_1))-(max(max(img_1))-min(min(img_1)))/2;
%img_1 = abs(img_1);
%img_1 = 1- img_1/max(max(img_1));
%vimg2 = vimg2/2 + img_1(1:kvad,1:kvad)/2;
%img_1n = (img_1-min(min(img_1)))/(max(max(img_1))-min(min(img_1)));
%img_1n(img_1n<=0.1)=1;
%vimg2 = vimg2./img_1n(1:kvad,1:kvad);
%tmp = 1 - abs(cos(pi*Gdir(1:kvad,1:kvad)/360));
%tmp = (Gdir(1:kvad,1:kvad)/140).^2;
%tmp = tmp./max(max(tmp))'
%vimg2 = vimg2/2 + tmp/2;
%vimg2 = vimg2./max(max(vimg2)); 
% img_1 = imgaussfilt(img_0,9);      
% tmp = img_1(1:kvad,1:kvad)./max(max(img_1(1:kvad,1:kvad)));
% vimg2 = vimg2.*(1-tmp);
%vimg2(vimg2 < 0.1) = eps;
% maxdx=max(max(ximg));
% maxdy=max(max(yimg));
% vimg=sqrt(ximg*ximg+yimg*yimg);
%vimg=ximg+yimg;
%maxdv=max(max(vimg));

%% Tikslo f-jos bazë - gradiento amplitudë
vimg2 = Gmag(1:kvad, 1:kvad);
vimg2 = vimg2./max(max(vimg2));                        % normavimas iki kvad
%figure; surf(vimg2);

%% STUFF
% xx = [floor(kvad/2):-1:0 1:floor(kvad/2)]'*ones(1,kvad);
% yy = xx';
% xy = sqrt(xx.^2+yy.^2);
% xy = xy(:);

%% Tikslo f-ja - L0 norma nuo gradientø skirtumo
vimg2 = 1000*vimg2(:);
for i = 1:kvad2
  for j = 1:kvad2;
    if(i==j)
      d(i,j) = 10;
    else
      d(i,j) = abs((vimg2(i)-vimg2(j)));%d(i,j) = sqrt((xy(i)-xy(j))^2 + (vimg2(i)-vimg2(j))^2);
    end
  end
end

%% Tikslo f-jos inversija & ...
vimg2 = max(d(:)) - d;
% ... normavimas &
vimg2 = vimg2/max(vimg2(:));
% ... slenkstinis filtravimas
vimg2(vimg2<0.85) = 0;

%% MHSA pradþia
c = vimg2;                  % - weight of city pair (C_i, C_j)
                            % t = 0, c_ij(t = 0) = max(d_ij) ? d_ij
                            % t > 0, c_ij in [0,1]
                            % larger c_ij , the stronger the mosquito m_ij
                            % and the easier it seeks and attacks the host
c(c==0) = eps;
r = 0.2*ones(size(c));      % - distance between mosquito m_ij and the host
                            % r_ij = 1 represents that the artificial mosquito m_ij is attacking the host and
                            % the shortest path passes through this host.
x = ones(size(c));          % - mosquito sex: 0 for male and 1 for female

%%
for i=1:kvad2
    c(i,:)=2*c(i,:)./(sum(c(i,:)));
end
for i=1:kvad2
    c(:,i)=2*c(:,i)./(sum(c(:,i)));
end
%d=[10 20 30 40 50 60 70 80 90 100];
%Pilno kelio skaiciavimas
%Z = 0;
%for i=1:kvad2
%    for j=1:kvad2
%        %d(atstumas)*r(praejimo koeficientas)
%        Z = Z + vimg2(i,j)*r(i,j);
%    end
%end
%s=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
%ts=[2 3 4 5 6 7 8 9 10 1 16 17 18 19 20 1 1 1 1 1];
%%

t = 1;
tic
while t < 1000,
  expuepclon = 0;
  sumrx = 0;
  rx = r.*x;
  %kiek kartu miestas sujungtas su kitais miestais(TSP problemoje 2 jungimai)
  for i=1:kvad2
    sumrx = sumrx + (sum(rx(i,:))-2);
  end
  for i=1:kvad2
    for j=1:kvad2
      %Pagal literatura u atstumas(radial distance) tarp miestu ir moskitu
      u(i,j) = 1-exp(-c(i,j)*r(i,j)*x(i,j));
      expuepclon = expuepclon+exp(-((u(i,j)*u(i,j))/(2*epclon*epclon)));
      %Pagal lit moskito judesys pagal
      dur(i,j) = -c(i,j)*x(i,j)*exp(-c(i,j)*r(i,j)*x(i,j));
      duc(i,j) = -r(i,j)*x(i,j)*exp(-c(i,j)*r(i,j)*x(i,j));
      dJr(i,j) = -c(i,j)*x(i,j)*exp(-c(i,j)*r(i,j)*x(i,j));
      dJc(i,j) = -r(i,j)*x(i,j)*exp(-c(i,j)*r(i,j)*x(i,j));
      dQr(i,j) = 2*sum(x(i,:))*sumrx+((1/(1+exp(-10*u(i,j))))-0.5)*dur(i,j);
      dQc(i,j) = ((1/(1+exp(-10*u(i,j))))-0.5)*duc(i,j);
    end
  end
  
  for i=1:kvad2
    for j=1:kvad2
      %Traukos funkcija prie uzkrestojo(host)
      dPr(i,j) = u(i,j)*(exp(-(u(i,j)*u(i,j))/(2*epclon*epclon))/expuepclon)*dur(i,j);
      dPc(i,j) = u(i,j)*(exp(-(u(i,j)*u(i,j))/(2*epclon*epclon))/expuepclon)*duc(i,j);
    end
  end
  %7-8 formules
  for i=1:kvad2
    for j=1:kvad2
      detr(i,j) = -lemta1*dur(i,j)-lemta2*dJr(i,j)-lemta3*dPr(i,j)-lemta4*dQr(i,j);
      detc(i,j) = -lemta1*duc(i,j)-lemta2*dJc(i,j)-lemta3*dPc(i,j)-lemta4*dQc(i,j);
    end
  end
  r = r + detr;
  c = c + detc;
  %gesinimas reiksmiu pagal blogiausia rezultata
  minr = min(min(r));
  minc = min(min(c));
  for i=1:kvad2
    for j=1:kvad2
      r(i,j) = r(i,j) - minr;
      c(i,j) = c(i,j) - minc;
    end
  end
  for i=1:kvad2
    r(i,:) = 2*r(i,:)./(sum(r(i,:)));
    c(i,:) = 2*c(i,:)./(sum(c(i,:)));
  end
  for i=1:kvad2
    r(:,i) = 2*r(:,i)./(sum(r(:,i)));
    c(:,i) = 2*c(:,i)./(sum(c(:,i)));
  end
  t = t + 1;
end
  %r1 = ones(size(r));
  %r1(r>0.8) = 0;
  
  img_01 = img_0(1:kvad,1:kvad);
  mx = (max(tril(r)));
  T = func_seperate_two_class(mx);
  img_01(reshape(mx,kvad,kvad)>T)=0;
  %img_01(r1==0)=0;
  %figure(101); surf(r); drawnow
  %disp(t)
  %figure(100); stem(max(tril(r))); %drawnow
  %figure(101); surf(img_01); view([ 0 90]); %drawnow ; pause(0.001)%view(-20,75); %pause

laikas=toc
%figure(101); surf(r)
%figure(102); surf(img_0)

%r1 = ones(size(r));
%r1(r>0.8) = 0;

%img_01 = img_0(1:kvad,1:kvad);
%img_01(r1==0)=0;

%figure(103); surf(img_01)

%vimg3 = r(1:kvad,1:kvad)
%[eps*ones(kvad,kvad) vimg2;...
%         vimg2 eps*ones(kvad,kvad)];
almost=ones(kvad).*255;
almost=almost-img_01;
ans = imresize(uint8(almost>=64).*255,4);
%disp(r);
%disp(' ');
%disp(c);



%Z = 0;
%figure(2)
%atv = img(10:30,30:50)
%Rezultatas(slenkstis pagal aki)
%for i=1:kvad
%    for j=1:kvad
%        pav(i,j)=uint8(img(i,j))*uint8((r(i,j)>0.2));
%    end
%end
%imshowpair(atv,pav,'montage');
