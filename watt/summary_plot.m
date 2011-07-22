load -ascii beam1.exx.dat;
nnn = rows(beam1_exx);
load -ascii beam1.xcoord.dat;
nx = rows(beam1_xcoord);
x = beam1_xcoord;
nt = nnn/nx;
for g=1:nt
    for h=1:nx
        ex(h,g) = beam1_exx((g-1)*nx+h,3);
    end
    time(g) = beam1_exx((g-1)*nx+1,2);
end
figure(1)
mesh (time,x,ex)
xlabel('time (s)')
ylabel('x (m)')
zlabel('Electric field (V/m)')

for g=1:nt
    efft(1:nx,g) = abs(fft(ex(1:nx,g)));
end

figure(2)
semilogy (time,efft(30,:),time,efft(40,:),time,efft(50,:))
xlabel('time (s)')
ylabel('Electric field amplitude (V/m)')


load -ascii beam1.fe0.dat;
mmm = rows(beam1_fe0);
load -ascii beam1.evcoord.dat;
nv = rows(beam1_evcoord);
vx = beam1_evcoord;
for g=1:nt
    for h=1:nv
        fe0(h,g) = beam1_fe0((g-1)*nv+h,3);
    end
end
figure(3)
plot(vx,fe0(:,1),vx,fe0(:,nt))
axis ([-5.0e+06,5.0e+06,0.0,0.8])
xlabel('v (m/s)')
ylabel('Distribution function (m^{-6}s^3)')
