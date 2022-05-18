% =================== %
% ISS Projekt 2021/22 %
% =================== %

close all; clear; clc;

%% 4.1 Základy
[audiodata, vzor_frekv] = audioread(['https://www.fit.vutbr.cz/' ...
    'study/courses/ISS/public/proj2021-22/signals/xvalen29.wav']);

delka = length(audiodata);
delka_sek = length(audiodata) / vzor_frekv;
MAX = max(abs(audiodata));
MIN = min(abs(audiodata));
mid = mean(audiodata);

disp("Vzorkovací frekvence = " + vzor_frekv + " Hz"); 
disp("Délka signálu ve vzorcích = " + delka);
disp("Délka signálu v sekundách = " + delka_sek + " s");
disp("Minimální hodnota normalizovaného signálu = " + MIN);
disp("Maximální hodnota normalizovaného signálu = " + MAX);

% Původní signál
delka_vzorek = delka_sek / (delka - 1);
t = 0 : delka_vzorek : delka_sek;

% Zobrazení signálu
figure;
plot(t, audiodata);
xlim([0 delka_sek]);
xlabel('Čas [s]');
ylabel('Amplituda');
title('Původní normalizovaný signál');

%% 4.2 Předzpracování a rámce
pocet_ramcu = floor(delka / 512) - 1;
vzorku_ramec = 1024;
prekryti = 512;
znely_ramec = 11;

% Rozdělení do rámců
ramec = zeros(vzorku_ramec, pocet_ramcu);
for index = 0 : pocet_ramcu - 1
    ramec(:, index + 1) = reshape(audiodata( ...
        index * prekryti + 1 : index * prekryti + vzorku_ramec ...
        ), vzorku_ramec, []);
end

% Znělý rámec
delka_ramce_sek = delka_sek / pocet_ramcu * 2;
t = 0 : delka_ramce_sek / (vzorku_ramec - 1) : delka_ramce_sek;

% Zobrazení znělého rámce
figure;
plot(t, ramec(:, znely_ramec));
xlim([0 delka_ramce_sek]);
xlabel('Čas [s]');
ylabel('Amplituda');
title('Znělý rámec');

%% 4.3 DFT
% Vytvoření W matice
omega = exp(-2 * pi * 1i / vzorku_ramec);
n = 0 : vzorku_ramec - 1;
[xx, yy] = meshgrid(n);
W = omega .^ (xx .* yy);

% DFT
DFT_vzorek = vzor_frekv / 2 / (vzorku_ramec / 2 - 1);
DFT = abs(W * ramec(:, znely_ramec));
DFT = DFT(1 : vzorku_ramec / 2);
t = 0 : DFT_vzorek : vzor_frekv / 2;

% Vestavěná FFT funkce
FFT = abs(fft(ramec(:, znely_ramec)));
FFT = FFT(1 : vzorku_ramec / 2);
frekvence = 0 : DFT_vzorek : vzor_frekv / 2;

% Porovnání DFT a FFT graficky
figure;
subplot(2, 1, 1);
plot(t, DFT);
xlim([0 vzor_frekv / 2]);
xlabel('Frekvence [hz]');
title('DFT funkce');

subplot(2, 1, 2);
plot(frekvence, FFT);
xlim([0 vzor_frekv / 2]);
xlabel('Frekvence [hz]');
title('Vestvěná FFT funkce');

% Porovnání transformací číselně
disp("Průměrný rodzíl hodnot = " + mean(FFT - DFT));
disp("Celkový rodzíl hodnot = " + sum(abs(FFT - DFT)));

%% 4.4 Spektogram
figure;
spectrogram(audiodata, vzorku_ramec, prekryti, [], vzor_frekv, 'yaxis'); 
colorbar;

xlabel('Čas [s]');
ylabel('Frekvence [kHz]');
title('Spektrogram');

%% 4.5 Určení rušivých frekvencí
% FFT pro první rámec
FFT = abs(fft(ramec(:, 1)));
FFT = FFT(1 : vzorku_ramec / 2);

figure;
plot(frekvence, FFT);
xlim([0 vzor_frekv / 2]);

xlabel('Frekvence [hz]');
title('FFT pro první rámec');

% Určení největších hodnot a jejich indexů
n = 1 : vzorku_ramec / 2;
FFT = [transpose(FFT); n];
sorted = sortrows(FFT.', 1).';
max_index = sort(sorted(2, end -3 : end));
max_index = frekvence(max_index);

% Rušivé frekvence
format ShortG; 
format compact;
n = 1 : 4;

disp("Rušivé frekvence = ");
disp(max_index);
max_index = max_index(1) * n;
disp(max_index);

%% 4.6 Generování signálu
linear = linspace(0, delka_sek - delka_vzorek, vzor_frekv * delka_sek + 1);
linear(end) = [];
linear = sin(transpose(max_index * 2 * pi) * linear);
linear = sum(linear) / 4;
audiowrite("audio/4cos.wav", linear, vzor_frekv);

% Porovnání spektrogramů
figure;
subplot(2, 1, 1);
spectrogram(audiodata, vzorku_ramec, prekryti, [], vzor_frekv, 'yaxis'); 
colorbar;
xlabel('Čas [s]');  
ylabel('Frekvence [kHz]');
title('Spektrogram původního signálu');

ax = gca;

subplot(2, 1, 2); 
spectrogram(linear, vzorku_ramec, prekryti, [], vzor_frekv, 'yaxis'); 
colorbar;
caxis([ax.CLim]);
xlabel('Čas [s]');
ylabel('Frekvence [kHz]');
title('Spektrogram rušivých frekvencí');

% Porovnání zvuků
sound(audiodata, vzor_frekv);
pause(delka_sek);
sound(linear, vzor_frekv); 
pause(delka_sek);

%% 4.7 Čistící filtr
% Hodnoty
zaverna_sirka = 30 / 2;
propustna_sirka = 50 + zaverna_sirka;
zvlneni = 3;
potlaceni = 40;
nyquist = vzor_frekv / 2;

% Vytvoření filtrů
B = zeros(4, 9);
A = zeros(4, 9);
H = dfilt.df2t.empty(0, 4);
figure

for index = 1 : 4
    x = max_index(index);

    % Bandstop butter filtr
    h = fdesign.bandstop(x - propustna_sirka, x - zaverna_sirka, ...
        x + zaverna_sirka, x + propustna_sirka, zvlneni, potlaceni, ...
        zvlneni, vzor_frekv);
    Hd = design(h, 'butter', 'MatchExactly', 'stopband');
    
    % Převedení hodnot do formátu [b, a]
    sos = get(Hd,'sosMatrix');
    [b, a] = sos2tf(sos);

    % Vypsání koeficientů
    format ShortG; 
    format compact
    disp("Koeficienty " + index + ":");
    disp(b);
    disp(a);

    % Uložení hodnot
    B(index, :) = b;
    A(index, :) = a;

    % Vykreslení impulzní odezvy
    subplot(2, 2, index);
    impz(b, a, 64, vzor_frekv);

    % Převedení do objektu
    H(index) = dfilt.df2t(B(index, 1 : end), A(index, 1 : end));
end

%% 4.8 Nulové body a póly
figure;
for index = 1 : 4
    subplot(2, 2, index);
    b(1 : end) = B(index, 1 : end);
    a(1 : end) = A(index, 1 : end);
    zplane(b, a);
end

%% 4.9 Frekvenční charakteristika
% Spojení filtrů do kaskády
Hcas = dfilt.cascade(H);
freqz(Hcas, 2^18, vzor_frekv);

%% 4.10 Filtrace
signal_filt = filter(Hcas, audiodata);

% Uložení audia
linear = linspace(0, delka_sek - delka_vzorek, vzor_frekv * delka_sek + 1);
linear(end) = [];
linear = sin(transpose(max_index * 2 * pi) * linear);
linear = sum(linear) / 4;
audiowrite("audio/clean_bandstop.wav", signal_filt, vzor_frekv);

%% Závěrečná srovnání
% Porovnání audia
sound(audiodata, vzor_frekv);
pause(delka_sek + 1);
sound(signal_filt, vzor_frekv);
pause(delka_sek);

% Porovnání spektrogramů
figure;
subplot(2, 1, 1);
spectrogram(audiodata, vzorku_ramec, prekryti, [], vzor_frekv, 'yaxis'); 
colorbar;
xlabel('Čas [s]');  
ylabel('Frekvence [kHz]');
title('Spektrogram původního signálu');

ax = gca;

subplot(2, 1, 2); 
spectrogram(signal_filt, vzorku_ramec, prekryti, [], vzor_frekv, 'yaxis'); 
colorbar;
caxis([ax.CLim]);
xlabel('Čas [s]');
ylabel('Frekvence [kHz]');
title('Spektrogram bez rušivých frekvencí');

% FFT pro poslední rámec
FFT = abs(fft(signal_filt(vzor_frekv-1024 : vzor_frekv)));
FFT = FFT(1 : vzorku_ramec / 2);

figure;
plot(frekvence, FFT);
xlim([0 vzor_frekv / 2]);

xlabel('Frekvence [hz]');
title('FFT pro poslední rámec bez rušivých frekvencí');
