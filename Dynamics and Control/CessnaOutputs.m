out = sim('Cessna');

t = out.simout.Time;
y = out.simout.Data;

y = squeeze(y);

if size(y,1) ~= length(t)
    y = y.';
end

figure
plot(t,y,'LineWidth',1.2)
grid on
xlabel('Time')
ylabel('Outputs')
legend('V','alpha (deg)','gamma (deg)','h','p','q')