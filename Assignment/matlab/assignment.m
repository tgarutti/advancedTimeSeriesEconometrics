%% Clear
clear
pad = '/Users/Rutger-Jan/Desktop/Assignment/Matlab' 
%pad = 'C:/Users/61849rla/Desktop/Assignment/Matlab';
% Change the above line to wherever you have saved this folder
cd(pad)

%% Load data

load('data')
dates = data(1:end,1)+datenum('31/12/1899','dd/mm/yyyy'); 
% column 1 contains date numbers
% we have to make a correction to the dates because Excel dates start counting at 1 Jan 1900
y = data(1:end,[2:3])'; % columns 2,3 and 4 contain the real data
% make sure y(t) is a column vector, so t runs towards the right

%% Plot series 1: GDP growth
plot(dates,y(1,:),'k')
dateFormat = 'yy';
datetick('x',dateFormat)
startrow=datenum('1/4/1947','dd/mm/yyyy')
endrow=datenum('1/4/2018','dd/mm/yyyy')
axis([startrow,endrow,-inf,inf])
title('Annualised GDP growth','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')

%% Plot series 2: Annualised inflation
plot(dates,y(2,:),'k')
dateFormat = 'yy';
datetick('x',dateFormat)
startrow=datenum('1/4/1947','dd/mm/yyyy')
endrow=datenum('1/4/2018','dd/mm/yyyy')
axis([startrow,endrow,-inf,inf])
title('Annualised inflation','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
