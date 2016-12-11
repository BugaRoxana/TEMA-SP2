%Formula de calcul a coeficientiilor seriei Fourier complexa au fost
%calculati manual.
%In conformitate cu cursul de SS un semnal se poate descompune intr-o suma
%de sinusuri si cosinusuri. Aces fapt este realizat prin intermediul
%seriilor Fourier. Ecuatia ce decrie o ramura ascendenta sau descendenta a
%unui semnal triunghiular poate fi usor descris cu ajutorul functiilor
%sinus si cosinus pe cand varfurile semnalului triunghiular ( inclusiv
%colturile de la nivelul de semnal 0) necesita un numar mai mare de
%componente spectrale intrucat un sinus sau un cosinus devine atat de
%ascutit doar la frecvente foarte inalte, ce sunt corespunzatoare unor
%componente spectrale de ordin ridicat. Pentru un semnal dreptunghiular,
%datorita pantelor extrem de abrupte, semnalele sinusoidate ce definesc
%aceste pante nu au timp sa se redreseze pentru zona pe palier aparand
%astfel niste spik-uri.

P=40;

D=4;

t=-120:0.1:120;

s=zeros(1,length(t));

m=1/4;

j=0;

for i=40:0.1:160
    j=j+1;
    c=floor(i/P);
    i=i-c*P;
    
    if(j==1)
        s(1,j+1200)=m*(t(j+1200));
    else if (i<=4)
            s(1,j+1200)=m*(t(j+1200))-m*40*(c-1);
            %Am adunca 40*(c-1) pentru a modela si celelalte perioade
        else if (i>=4&&i<=8)
                s(1,j+1200)=-m*(t(j+1200))+2+m*40*(c-1);
            end
        end
    end
end

j=0;

for i=40:0.1:160
    j=j+1;
    c=floor(i/P);
    i=i-c*P;
    
    if(j==1)
        s(1,j)=m*(t(j+1200));
    else if (i<=4)
            s(1,j)=m*(t(j+1200))-m*40*(c-1);
            %Am adunca 40*(c-1) pentru a modela si celelalte perioade
        else if (i>=4&&i<=8)
                s(1,j)=-m*(t(j+1200))+2+m*40*(c-1);
            end
        end
    end
end


N=50;

p=-9:1:9;

Ak=zeros(length(2*p));
Ak(10)=0.1;

for k=1:2:9
    
        %Ak(k+10)=4*0.5/(k*k*pi*pi);
        %Ak(10-k)=Ak(k+10);
        %Formula calculata dupa calculul lui Ck
        Ak(k+10)=5/(k*k*pi*pi)*abs((2*(cos(k*pi/5)-1j*sin(k*pi/5))-(cos(k*2*pi/5)-1j*sin(k*2*pi/5))-1));
        Ak(10-k)=Ak(k+10);
end

w0=2*pi/P;

r=linspace(-120,120,241);
f=0*r;

for k=-N:1:N
    
   if (k==0)
       continue;
   end
   
   C_k=2.5/(k*k*pi*pi)*(2*(cos(k*pi/5)-1j*sin(k*pi/5))-(cos(k*2*pi/5)-1j*sin(k*2*pi/5))-1);
   f_k=C_k*exp(1j*k*w0*r);
   f=f+f_k;
   
    
end

f=f+0.1;

figure (1)

plot(t,s),grid,xlabel('Timp'),ylabel('s[t]'),title('Graficul lui s in functie de timp');

figure (2)

stem(p,Ak),grid,xlabel('k'),ylabel('Ak'),title('Spectrul de amplitudini ale lui s');

figure (3)

plot(r,f),grid,xlabel('f[t]'),ylabel('Timp'),title('Graficul dependentei functiei s de timp refacut cu ajutorul seriei Fourier exponentiala');

