clear all
clc
% % Motor de Corriente Continua con Carga % %

% Parametros Tecnicos del Motor
Vnom=9;
R=0.113;
L= 0.064e-3;
Bm=30/(pi*9.19*1000);
Kb = 1/(881*2*pi/60);
Ki=10.8e-3;
Jm = 12.8e-7;
Gain=5;

%Acondicionador de Senales
VoltajeOtorgado=10;
Ga=VoltajeOtorgado/(2);


%Carga%

%Ruedas
Rrueda=(0.03); %diametro de 10 cm
MasaTotal=4; %Masa de ruedas mas la masa del auto y del motor. El original era 16kg.
Jc=(MasaTotal/2)*(Rrueda^2); %momento de inercia de la carga
Bc=0;


noloadspeed=6430;
stalltorque=0.7;
%Engranajes
MAXVelocidadAuto=40; %km/hora
pendiente = -stalltorque/noloadspeed;
gravedad=9.81;
rpm=(MasaTotal*gravedad*0.005-stalltorque)/pendiente;
wmotor=rpm*2*pi/60;
wauto=MAXVelocidadAuto/3.6/Rrueda;
RelacionEngranajes = wauto/wmotor;
RelacionEngranajes = 15/20; %Aproximacion

%Ejes
Reje = 0.005; %diametro de un centimetro
Meje=0.02; %Pesa 20 gramos
Je = Meje * (Reje^2); %momento de inercia de los ejes

% Momento de Inercia y Rozamiento Viscoso con Carga
Jt=Jm+Jc*(RelacionEngranajes^2); %Momento de inercia total
Bt=Bm+Bc;

%Problematica (en metros)
PosInicial=0;
Objetivo=1;

% Especificaciones Transistorio
tr=1.4; % Tiempo de Levantamiento (lo que tarda en llegar del 10% al 90%). El original era 10
Mp=0.1/100; %El sobre paso no puede ser 0! El original era 2/100


zita=1/sqrt(((-pi/log(Mp))^2)+1);
%if(zita<0.69)
%    wn=3.2/(zita*ts); 
%else
%    wn=4.5*zita/ts; 
%end
    
    
wn=(0.8+2.5*zita)/tr; 


F_Ideal = tf(wn^2, [1, 2*zita*wn, wn^2]);
display('Sistema "Ideal"')
stepinfo(F_Ideal)
polosDeseados = [-zita*wn+1i*wn*sqrt(1-zita^2);-zita*wn-1i*wn*sqrt(1-zita^2)];
p1 = polosDeseados(1);
p2 = polosDeseados(2);


%-------------Sistema Sin Compensar-------------%

I = tf([1],[1 0]); %Integrador
Gmotor = tf([Ki],[L*Jt R*Jt+Bm*L R*Bm+Ki*Kb]); %FT del motor
Gaux = Gmotor * I * RelacionEngranajes*Rrueda; %G
Sensor = tf(10/2, [0.1 1]);%Encoder [0.1 1]

FTLA = Gaux*Sensor;
FTLC = Ga*feedback(Gaux, Sensor);
%-------------Bisectriz-------------%
[z,p,k]=zpkdata(FTLA,'v');
%aporte de fase de los polos al punto de diseño
angp1=180-atand(abs(imag(p(1))-imag(p1))/abs(real(p1)-real(p(1))));    
angp2=atand(abs(imag(p(2))-imag(p1))/abs(real(p1)-real(p(2))));
angp3=atand(abs(imag(p(3))-imag(p1))/abs(real(p1)-real(p(3))));
angp4=atand(abs(imag(p(4))-imag(p1))/abs(real(p1)-real(p(4))));

%angulo del compensador
anguloAportado=-angp1-angp2-angp3-angp4;
angComp=-180-anguloAportado;
angComp=abs(angComp);

%compensador
angbisectriz = 180-acosd(zita);
augaux1=180-(angbisectriz/2-angComp/2)-acosd(zita);
augaux2=180-(angbisectriz/2+angComp/2)-acosd(zita);
cero=-sind(angbisectriz/2-angComp/2)*abs(p1)/sind(augaux1);
polo=-sind(angbisectriz/2+angComp/2)*abs(p1)/sind(augaux2);
C=zpk(cero,polo,1);

%ganacia del compensador
modp1=abs(p1-p(1));
modp2=abs(p1-p(2));
modp3=abs(p1-p(3));
modp4=abs(p1-p(4));
modp=abs(p1-polo);
modz=abs(p1-cero);
K=(modp1*modp2*modp3*modp4*modp)/(modz*k);
Comp = K*C;

%-------------respuesta en frecuencia-------------
FTLC_Compensada = Ga*feedback(Gaux*Comp, Sensor);
display('Sistema sin Compensar')
stepinfo(FTLC)



%Tabla Routh-Hurwitz
EcuacionCaracteristica = 1+Gaux*Sensor;
x = sym ('k');
A=6.488192000000000e-09;
B=1.152724613742626e-05;
C=1.380717209807807e-04;
D=2.344807880651816e-04;
E=0.001215000000000*x;
RH(1).Col1=A;
RH(1).Col2=C;
RH(1).Col3=E;
RH(2).Col1=B;
RH(2).Col2= D;
RH(2).Col3=0;
RH(3).Col1= -(A*D-B*C)/B;
RH(3).Col2=-(A*0-E*B)/B;
RH(3).Col3=0;
RH(4).Col1=-(B*(RH(3).Col2)-D*(RH(3).Col1))/(RH(3).Col1); %*k
RH(5).Col1=RH(3).Col2;
LimiteSuperior=2.309376963236687;
FCritica=Ga*feedback(Gaux*LimiteSuperior, Sensor);
%Imprimir Tabla
lalala=1;
while(lalala<4)
    RH(lalala).Col1;
    RH(lalala).Col2;
    lalala=lalala+1;
end

display('Sistema Compensado por el Metodo de Lugar de Raices (con ganancia ajustada)')
stepinfo(FTLC_Compensada)

[z_ftlc,p_ftlc,k_ftlc]=zpkdata(FTLC,'v');
Ess_compensado = FTLA*Comp;

%Por alguna puta razon esta dando mal el error en regimen permanente en
%Simulink!!!!
%Compensador en Atraso
KvCompensado=1/0.15;
KvSinCompensar=0.5819;
Beta = KvCompensado/KvSinCompensar;
PoloAtraso=0.01;
T=1/(PoloAtraso*Beta);
CeroAtraso=1/T;
CompensadorAtraso = zpk(-CeroAtraso, -PoloAtraso, 1);
FTLC_Super_Compensada = Ga*feedback(CompensadorAtraso*Comp*Gaux, Sensor);

pole(FTLC_Compensada)
polosDeseados
stepinfo(FTLC_Super_Compensada)

%Calculo del error verdadero
Gcomp = Comp*Gaux;
Geq=Gcomp*Ga/(1+Gcomp*(Sensor-Ga));
%Compensacion del error verdadero
%KvSinCompensar = 19.9997;
%KvCompensado=100;
%Beta = KvCompensado/KvSinCompensar;
%PoloAtraso2=0.01;
%T2=1/(PoloAtraso*Beta);
%CeroAtraso=1/T;
%CompensadorAtraso = zpk(-CeroAtraso, -PoloAtraso, 1);
