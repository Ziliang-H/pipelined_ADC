%8位1-bit/stage Pipelined ADC建模
clear; clc; close all;

%参数设置
N_total = 8; 
V_ref = 1; 
Fs = 10e6; 
Fin = 0.123e6; 
N_sample = 2^14;
Vin_amp = 0.45*V_ref;
t = (0:N_sample-1)/Fs;
Vin = Vin_amp*sin(2*pi*Fin*t) + V_ref/2;

%8级流水线处理
bit_all = zeros(N_total, N_sample);
residue = Vin;
for stage = 1:N_total
    [bit_all(stage,:), residue] = stage_1bit(residue, V_ref);
end

%模拟输出
Vout = 0;
for stage = 1:N_total
    Vout = Vout + bit_all(stage,:) * 2^(N_total - stage);
end
Vout = Vout * V_ref / (2^N_total - 1);

%时域图
figure(1);
plot(t(1:100), Vin(1:100), 'b-', t(1:100), Vout(1:100), 'r.-');
xlabel('时间(s)'); ylabel('电压(V)'); legend('输入','输出');
title('8位流水线ADC输入输出时域'); grid on;

%频谱图
figure(2);
n_half = floor(N_sample/2)+1;
Vin_fft = abs(fft(Vin))/N_sample; 
Vin_fft = Vin_fft(1:n_half);
Vout_fft = abs(fft(Vout))/N_sample; 
Vout_fft = Vout_fft(1:n_half);
Vin_fft(2:end) = 2*Vin_fft(2:end); 
Vout_fft(2:end) = 2*Vout_fft(2:end);
f_pos = linspace(0, Fs/2, n_half);
plot(f_pos, 20*log10(Vout_fft/max(Vin_fft)), 'r-', ...
    f_pos, 20*log10(Vin_fft/max(Vin_fft)), 'b-');
xlabel('频率(Hz)'); 
ylabel('幅度(dB)'); 
legend('输出','输入');
title('ADC输入/输出频谱'); 
xlim([0 Fs/4]); 
ylim([-80 10]); 
grid on;

%residue
figure(3);
Vin_test = linspace(0, V_ref, 1000);
[~, ~, G_actual, Vth_actual, Vref_coeff] = stage_1bit(0, V_ref); % 调用函数获取参数
%实际
res_act = arrayfun(@(v)stage_residue(v, V_ref, G_actual, Vth_actual, Vref_coeff), Vin_test);
%理想
Vth_ideal = V_ref/2; 
G_ideal = 2;
res_ideal = (Vin_test<Vth_ideal).* (G_ideal*Vin_test) + ...
            (Vin_test>Vth_ideal).* (G_ideal*Vin_test - V_ref);
% 绘图
plot(Vin_test, res_act, 'b-', 'LineWidth', 1.5); 
hold on;
plot(Vin_test, res_ideal, 'k--', 'LineWidth', 1.2);
xline(Vth_actual, 'r--', '实际阈值'); 
xline(Vth_ideal, 'g:', '理想阈值');
yline(V_ref, 'm-.', '理论峰值'); 
yline(0, 'c-.', '理论谷值');
xlabel('Vin (V)'); ylabel('Residue (V)'); 
title(sprintf('1-bit/stage余量函数 G=%.4f, Vth=%.4fV', G_actual, Vth_actual));
legend('实际','理想','Location','best'); 
grid on;
hold off;


function [bit_out, residue, G_actual, Vth_actual, Vref_coeff] = stage_1bit(Vin, V_ref)
% 电容mismatch
mismatch_cap = 0.05; 
C1 = 1 + mismatch_cap * randn();
C2 = 1 + mismatch_cap * randn();
G_cap = (C1 + C2)/C1;
% gain error
A_opamp = 1000; % 运放开环增益
A_actual = A_opamp * (1 + 0.05 * randn()); % 增益偏差
G_actual = G_cap / (1 + G_cap/A_actual);
Vref_coeff = (C2/C1) / (1 + 1/A_actual);
% 比较器失调
V_os = V_ref * 0.05 * randn(); 
Vth_actual = V_ref/2 + V_os;

bit_out = double(Vin > Vth_actual);
residue = G_actual * Vin - bit_out * Vref_coeff * V_ref;
end

%根据给定参数计算单点余量
function r = stage_residue(v, V_ref, G, Vth, Vref_coeff)
if v > Vth
   r = G * v - Vref_coeff * V_ref;
else
   r = G * v;
end
end

