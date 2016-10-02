function [xp] = Buckley()
% Tempo (d):
t = 0.2;
dt = 0.005;

%% ------------------------------- Geometria ------------------------------
% Comprimento (m):
L = 1;
% Área da secção transversal (m²):
A = 1;
%% ------------------ Propriedades dos fluidos e da rocha -----------------
% Porosidade (%):
phi = 1;
% Viscosidades (cp):
muo = 1;
muw = 1;
% Fator volume da formação da água (m³/m³ std):
Bw = 1;
%% ------------------- Condições iniciais e de contorno -------------------
% Saturação de óleo residual (%):
sor = 0.1;
% Saturação de água (%) (condição inicial):
swi = 0.1;
% Vazão de injeção de água (m³ std/d) (Como é std tem que usar o Bw):
Q = 1;

% Criando variáveis simbólicas (para poder mudar o modelos):
syms sw u dfwdsw fw krw kro swn;
% Permeabilidades relativas (Usando o modelo de Corey):
no = 2;
nw = 2;
% Saturação normalizada:
swn = (sw-swi)./(1-swi-sor); %#ok<NODEF> 
kro = (1-swn).^no;
krw = swn.^nw;

% Fluxo fracionário (lambda_w/lmabda_t):
fw = (krw/muw)/((krw/muw)+(kro/muo));
% derivada simbolica da função de fluxo:
dfwdsw = diff(fw,sw)

% Weldge tangent saturation (Le Veque 2002) - Resolve a equação para achar 
% a saturação da frente:
st = solve((sw-swi)*dfwdsw-fw,sw);
roots = double(st);
roots = roots(roots<=1);
roots = roots(roots>=0);
swf = max(roots);
% Saturação da frente:
sw = swf; %#ok<NASGU>

% Função saturação em função da posição:
%
%           / Sw(i) = Sw(x) = (inversa de x(sw)), se x(i) <= xf
%   Sw(i) = |
%           \ Sw(i) = swi, se x(i) >= xf

% Como ja podemos calcular o valor da posição que em determinado tempo terá
% uma determinada satuação criamos um vetor com espaçamento constante de 
% saturações, calculamos o x correspondente e plotamos (x,sw).
sw = (0:0.01:1)';

% Considerando-se que a vazão é constante no tempo e igual a Q a injeção
% acumulada é dada por: W = Q*t. Desta forma pode-se escrever a posição
% como um perfil (multiplo da derivada do fluxo) vezes o tempo
% Perfil:
% Equação para a posição (dividida pelo tempo), Eq. (14.214) do livro 
% Engenharia de reservatórios de petróleo (Livro verde do Adalberto):

%Descomente para ver os graficos da funções
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
    % Atualizando o perfil (se a vazão de injeção for variável é 
    % necessário mudar a formula)
    xp = profile*time;
    time = time + dt;
    grid on;
    hold on;
    % O perfil só tem valores com significado fisico quando sw >= swf, logo
    % plotamos apenas esses pontos:
    plot(xp(sw>=swf),sw(sw>=swf),'--b','LineWidth',1.5);
    % Abaixo do valor de swf o que se tem é o choque e a saturação inicial
    % sendo assim completamos o perfil da função no gráfico concatenando e
    % plotando esses valores:
    xf = max(xp(sw>=swf));
    plot([xf xf L],[swf swi swi],'--b','LineWidth',1.5);
    axis([0 L 0 1]);
    xlabel('Posição (m)');
    ylabel('Sw');
    title(['Perfil de saturações para o tempo (t = ' num2str(time) ' d)']);
    pause(0.1);
end

end