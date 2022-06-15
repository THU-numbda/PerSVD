clear;
fid = fopen('Dense1.dat', 'r');
A = fread(fid, [1000, 1000], 'float=>double');
A = A';
pmax = 4;

k=100;

[m, n] = size(A);
Omega = randn(n, 1.5*k);
[U, S, V] = svd(A, 'econ');
ss = diag(S(1:k, 1:k)).^2;
s101 = S(k+1, k+1).^2;
A_Ak = A - U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';
Ak_f = norm(A_Ak, 'fro');
Ak_2 = norm(A_Ak, 2);

err1_pve = [];
err1_f = [];
err1_s = [];
err2_pve = [];
err2_f = [];
err2_s = [];
err3_pve = [];
err3_f = [];
err3_s = [];
err4_pve = [];
err4_f = [];
err4_s = [];

for p = 0:1:pmax
    
    [u1, s1, v1] = basic_rSVD(A, k, p, Omega);
    sst = diag((u1'*A)*A'*u1);
    pvet = max(abs(sst-ss)./s101);
    err1_pve = [err1_pve, pvet];
    err1_f = [err1_f, (norm(A-u1*s1*v1', 'fro')-Ak_f)/Ak_f];
    err1_s = [err1_s, (norm(A-u1*s1*v1', 2)-Ak_2)/Ak_2];
    
    [u1, s1, v1] = rSVD_fp(A, k, p, Omega);
    sst = diag((u1'*A)*A'*u1);
    pvet = max(abs(sst-ss)./s101);
    err4_pve = [err4_pve, pvet];
    err4_f = [err4_f, (norm(A-u1*s1*v1', 'fro')-Ak_f)/Ak_f];
    err4_s = [err4_s, (norm(A-u1*s1*v1', 2)-Ak_2)/Ak_2];
    
    [u2, s2, v2] = PerSVD_once(A, k, p, Omega);
    sst = diag((u2'*A)*A'*u2);
    pvet = max(abs(sst-ss)./s101);
    err2_pve = [err2_pve, pvet];
    err2_f = [err2_f, (norm(A-u2*s2*v2', 'fro')-Ak_f)/Ak_f];
    err2_s = [err2_s, (norm(A-u2*s2*v2', 2)-Ak_2)/Ak_2];
    
    [u3, s3, v3] = PerSVD_update(A, k, p, Omega);
    sst = diag((u3'*A)*A'*u3);
    pvet = max(abs(sst-ss)./s101);
    err3_pve = [err3_pve, pvet];
    err3_f = [err3_f, (norm(A-u3*s3*v3', 'fro')-Ak_f)/Ak_f];
    err3_s = [err3_s, (norm(A-u3*s3*v3', 2)-Ak_2)/Ak_2];
end

x = 2:2:10;
figure(1);
semilogy(x, err1_pve, 's-',x./2, err4_pve, '+-', x./2, err2_pve, 'o-', x./2, err3_pve, '^-');
xlabel('#passes')
ylabel('\epsilon_{PVE}');
legend('Alg. 1','Alg. 3', 'Alg. 4', 'Alg. 5');
xmin = 1;
xmax = 10;
ymin = 0.5*min(err3_pve);
ymax = 2*err1_pve(1);
axis([xmin, xmax, ymin, ymax]);
figure_FontSize=25; 
set(findobj('FontSize',10),'FontSize',25); 
set( get(gca,'XLabel'),'FontSize',figure_FontSize); 
set( get(gca,'YLabel'),'FontSize',figure_FontSize); 
set( get(gca,'XAxis'),'FontSize',figure_FontSize); 
set( get(gca,'YAxis'),'FontSize',figure_FontSize); 
set( get(gca,'XAxis'),'LineWidth',2); 
set( get(gca,'YAxis'),'LineWidth',2); 
set( get(gca,'Legend'),'FontSize',20); 
set(findobj( get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca, 'YTick', [1e-8,1e-6, 1e-4, 1e-2, 1, 100]);

figure(2);
semilogy(x, err1_f, 's-', x./2, err4_f, '+-', x./2, err2_f, 'o-', x./2, err3_f, '^-');
xlabel('#passes')
ylabel('\epsilon_{F}');
legend('Alg. 1', 'Alg. 3', 'Alg. 4', 'Alg. 5');
xmin = 1;
xmax = 10;
ymin = 0.5*min(err3_f);
ymax = 2*err1_f(1);
axis([xmin, xmax, ymin, ymax]);
figure_FontSize=25; 
set(findobj('FontSize',10),'FontSize',25); 
set( get(gca,'XLabel'),'FontSize',figure_FontSize); 
set( get(gca,'YLabel'),'FontSize',figure_FontSize); 
set( get(gca,'XAxis'),'FontSize',figure_FontSize); 
set( get(gca,'YAxis'),'FontSize',figure_FontSize); 
set( get(gca,'XAxis'),'LineWidth',2); 
set( get(gca,'YAxis'),'LineWidth',2); 
set( get(gca,'Legend'),'FontSize',20); 
set(findobj( get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca, 'YTick', [1e-8,1e-6, 1e-4, 1e-2, 1, 100]);

figure(3);
semilogy(x, err1_s, 's-',x./2, err4_s, '+-', x./2, err2_s, 'o-', x./2, err3_s, '^-');
xlabel('#passes')
ylabel('\epsilon_{s}');
legend('Alg. 1', 'Alg. 3', 'Alg. 4', 'Alg. 5');
xmin = 1;
xmax = 10;
ymin = 0.5*min(err3_s);
ymax = 2*err1_s(1);
axis([xmin, xmax, ymin, ymax]);
figure_FontSize=25; 
set(findobj('FontSize',10),'FontSize',25); 
set( get(gca,'XLabel'),'FontSize',figure_FontSize); 
set( get(gca,'YLabel'),'FontSize',figure_FontSize); 
set( get(gca,'XAxis'),'FontSize',figure_FontSize); 
set( get(gca,'YAxis'),'FontSize',figure_FontSize); 
set( get(gca,'XAxis'),'LineWidth',2); 
set( get(gca,'YAxis'),'LineWidth',2); 
set( get(gca,'Legend'),'FontSize',20); 
set(findobj( get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca, 'YTick', [1e-8,1e-6, 1e-4, 1e-2, 1, 100]);