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
VoltajeOtorgado=0.5;
Ga=VoltajeOtorgado/(360);
Gb=3500/10*VoltajeOtorgado/360;

%Carga%

%Ruedas
Rrueda=(0.03); %diametro de 10 cm
MasaTotal=16; %Masa de ruedas mas la masa del auto y del motor.
Jc=(MasaTotal/2)*(Rrueda^2); %momento de inercia de la carga
Bc=0;


%Engranajes
MAXVelocidadAuto=40; %km/hora
wmotor=10000*2*pi/60;
wauto=MAXVelocidadAuto/3.6/Rrueda;
%RelacionEngranajes=1; %Esta relacion es de prueba.
RelacionEngranajes = wauto/wmotor; 

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
tr=1.4; % Tiempo de Levantamiento (lo que tarda en llegar del 10% al 90%)
Mp=0.1/100; %El sobre paso no puede ser 0!


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
DisAng = 360/(2*pi*Rrueda);
Potenciometro=10/(pi*3500/180);
I = tf([1],[1 0]); %Integrador
Gmotor = tf([Ki],[L*Jt R*Jt+Bm*L R*Bm+Ki*Kb]); %FT del motor
Gaux = Gmotor * I * RelacionEngranajes; %G
%SensorLaser=tf([10/5],[0.002 1]);
Sensor = Potenciometro * Gb; %H

FTLA = Gaux*Sensor;
FTLC = DisAng*Ga*feedback(Gaux, Sensor)*Rrueda;
%-------------Bisectriz-------------%
[z,p,k]=zpkdata(FTLA,'v');
%aporte de fase de los polos al punto de diseño
angp1=180-atand(abs(imag(p(1))-imag(p1))/abs(real(p1)-real(p(1))));    
angp2=atand(abs(imag(p(2))-imag(p1))/abs(real(p1)-real(p(2))));
angp3=atand(abs(imag(p(3))-imag(p1))/abs(real(p1)-real(p(3))));

%angulo del compensador
angComp=-180+angp1+angp2+angp3;

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
modp=abs(p1-polo);
modz=abs(p1-cero);
K=(modp1*modp2*modp3*modp)/(modz*k);
%K=1.3;
Comp = K*C;

%-------------Grafica de rlocus-------------%
%rlocus(FTOL);axis([-50000 3000 -10000 10000])
%rlocus(FTOL);axis([-5 5 -500 500])
%rlocus(FTOL*C);axis([-50000 3000 -10000 10000])
%rlocus(FTOL*C);axis([-5 5 -500 500]);
%FTOL2=Comp*FTOL;

%-------------respuesta en frecuencia-------------
FTLC_Compensada = DisAng*Ga*feedback(Gaux*Comp, Sensor)*Rrueda;
%step(FTLC)
display('Sistema sin Compensar')
stepinfo(FTLC)
%stepinfo(FTLC_Compensada)
%sys=FTCL2; 
%margin(FTCL2);grid on; 
%BW = bandwidth(FTCL2);
%nyquist(FTCL2);


%Tabla Routh-Hurwitz
EcuacionCaracteristica = 1+Gmotor * tf(1, [1, 0]) * RelacionEngranajes*Potenciometro*Gb;
x = sym ('x');
A=L*Jt;
B=R*Jt+Bm*L;
C=Ki*Kb+R*Bm;
RH(1).Col1=A;
RH(1).Col2=C;
RH(2).Col1=B;
RH(2).Col2=(Sensor*Ki*RelacionEngranajes); %*k
RH(3).Col1= ((B*C)-(A*Sensor*Ki*RelacionEngranajes*x))/B; %*k
RH(3).Col2=0;
RH(4).Col1=Sensor*Ki*RelacionEngranajes; %*k
LimiteSuperior=(B*C)/(A*Sensor*Ki*RelacionEngranajes);
FCritica=DisAng*Ga*feedback(Gaux*LimiteSuperior, Sensor)*Rrueda;
%Imprimir Tabla
lalala=1;
while(lalala<4)
    RH(lalala).Col1;
    RH(lalala).Col2;
    lalala=lalala+1;
end

display('Sistema Compensado por el Metodo de Lugar de Raices (con ganancia ajustada)')
stepinfo(FTLC_Compensada)

%%%%Compensacion por variables de estado%%%%%
display('Variables de estado')
[n, d] = tfdata(Gaux, 'v');     %Coeficientes del numerador y el denominador de la FTLC
[E1, E2, E3, E4] = tf2ss(n, d); %Obtengo el sistema en espacio de estados
%Determino la Observabilidad del sistema
Observabilidad = [E3; E3*E1; E3*(E1^2)];
if det(Observabilidad)~=0
    display('El sistema es observable.')
    Observable=true;
else
    display('El sistema no es observable.')
    Observable=false;
end
%Determino la Controlabilidad del sistema
Controlabilidad = [E2 E1*E2 (E1^2)*E2];
if det(Controlabilidad)~=0
    display('El sistema es controlable.');
    Controlable=true;
else
    display('Cuidado! El sistema no es controlable.')
    Controlable=false;
end

if(Controlable && Observable) %Si el sistema es controlable, sera compensado por variables de estado.
    alpha=5;  
    polosExtras= -(alpha*wn);  
                               %Agregamos un polo lo suficientemente
                               %alejados para que no interfieran con la
                               %respuesta del sistema.
    vectorK = place(E1, E2, [transpose(polosDeseados), polosExtras]); %Creo el vector k para el calculo del compensador
    Sistema_Compensado_Variables_Estado=ss(E1-E2*vectorK, E2, E3, E4); %Creo el sistema compensado en variables de estado
    [n_compensado, d_compensado] = ss2tf(E1-E2*vectorK, E2, E3, E4); %Lo transformo a TF
    FT_Compensado_Variables_Estado = tf(n_compensado, d_compensado);
    [z_estado, p_estado, k_estado] = tf2zpk(n_compensado, d_compensado);
    display('Sistema compensado por Variables de Estado')
    stepinfo(FT_Compensado_Variables_Estado)
end


%Criterio de polos dominantes, nos permite introducir otro polo de manera
%que este no afecte de manera significante a la dinamica del sistema.


