clear 
clc

%%


DSS = readmatrix("LDPC_SPA DSS.csv");
IRCMS = readmatrix("LDPC_SPA IRCMS.csv");
QC = readmatrix("LDPC_SPA QC_PEG.csv");




figure(1)

semilogy(IRCMS(1,:), IRCMS(2,:), 'b-pentagram');
hold on 
semilogy(DSS(1,:), DSS(2,:), 'r-pentagram');
hold on 
semilogy(QC(1,:), QC(2,:), 'k-pentagram');
hold on 

semilogy(IRCMS(1,:), IRCMS(3,:), 'b-o');
hold on 
semilogy(DSS(1,:), DSS(3,:), 'r-o');
hold on 
semilogy(QC(1,:), QC(3,:), 'k-o');
hold on 

legend('(1200,600) IRCMS:BER','(1200,600) DSS:BER','(1200,600) QC PEG:BER','(1200,600) IRCMS:FER','(1200,600) DSS:FER','(1200,600) QC PEG:FER');
xlabel('Eb/N0');
%ylabel('BER');
axis([1.5 3.5 10^-7 1]);

