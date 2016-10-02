function [xp] = Buckley()
% Tempo (d):
t = 0.2;
dt = 0.005;

%% ------------------------------- Geometria ------------------------------
% Comprimento (m):
L = 1;
% �rea da sec��o transversal (m�):
A = 1;
%% ------------------ Propriedades dos fluidos e da rocha -----------------
% Porosidade (%):
phi = 1;
% Viscosidades (cp):
muo = 1;
muw = 1;
% Fator volume da forma��o da �gua (m�/m� std):
Bw = 1;
%% ------------------- Condi��es iniciais e de contorno -------------------
% Satura��o de �leo residual (%):
sor = 0.1;
% Satura��o de �gua (%) (condi��o inicial):
swi = 0.1;
% Vaz�o de inje��o de �gua (m� std/d) (Como � std tem que usar o Bw):
Q = 1;

% Criando vari�veis simb�licas (para poder mudar o modelos):
syms sw u dfwdsw fw krw kro swn;
% Permeabilidades relativas (Usando o modelo de Corey):
no = 2;
nw = 2;
% Satura��o normalizada:
swn = (sw-swi)./(1-swi-sor); %#ok<NODEF> 
kro = (1-swn).^no;
krw = swn.^nw;

% Fluxo fracion�rio (lambda_w/lmabda_t):
fw = (krw/muw)/((krw/muw)+(kro/muo));
% derivada simbolica da fun��o de fluxo:
dfwdsw = diff(fw,sw)

% Weldge tangent saturation (Le Veque 2002) - Resolve a equa��o para achar 
% a satura��o da frente:
st = solve((sw-swi)*dfwdsw-fw,sw);
roots = double(st);
roots = roots(roots<=1);
roots = roots(roots>=0);
swf = max(roots);
% Satura��o da frente:
sw = swf; %#ok<NASGU>

% Fun��o satura��o em fun��o da posi��o:
%
%           / Sw(i) = Sw(x) = (inversa de x(sw)), se x(i) <= xf
%   Sw(i) = |
%           \ Sw(i) = swi, se x(i) >= xf

% Como ja podemos calcular o valor da posi��o que em determinado tempo ter�
% uma determinada satua��o criamos um vetor com espa�amento constante de 
% satura��es, calculamos o x correspondente e plotamos (x,sw).
sw = (0:0.01:1)';

% Considerando-se que a vaz�o � constante no tempo e igual a Q a inje��o
% acumulada � dada por: W = Q*t. Desta forma pode-se escrever a posi��o
% como um perfil (multiplo da derivada do fluxo) vezes o tempo
% Perfil:
% Equa��o para a posi��o (dividida pelo tempo), Eq. (14.214) do livro 
% Engenharia de reservat�rios de petr�leo (Livro verde do Adalberto):

%Descomente para ver os graficos da fun��es
% figure(1);
% grid on;
% hold on;
% plot(sw,double(eval(dfwdsw)),'.-b');
% plot(sw,double(eval(fw)),'.-r');
% axis([0 1 0 2]);
% hold off;


profile = (Bw*Q/(A*phi))*double(eval(dfwdsw));

%figure(2);
time = 0;
while (time<t)
    cla;
    % Atualizando o perfil (se a vaz�o de inje��o for vari�vel � 
    % necess�rio mudar a formula)
    xp = profile*time;
    time = time + dt;
    grid on;
    hold on;
    % O perfil s� tem valores com significado fisico quando sw >= swf, logo
    % plotamos apenas esses pontos:
    plot(xp(sw>=swf),sw(sw>=swf),'--b','LineWidth',1.5);
    % Abaixo do valor de swf o que se tem � o choque e a satura��o inicial
    % sendo assim completamos o perfil da fun��o no gr�fico concatenando e
    % plotando esses valores:
    xf = max(xp(sw>=swf));
    plot([xf xf L],[swf swi swi],'--b','LineWidth',1.5);
    axis([0 L 0 1]);
    xlabel('Posi��o (m)');
    ylabel('Sw');
    title(['Perfil de satura��es para o tempo (t = ' num2str(time) ' d)']);
    pause(0.1);
end

end