function plotGratingReflectanceSurf(lambdas, R_original, R_received, length_domain, dane_po_przetworzeniu, dane_w_kanale, dane_referencyjne, N_s,D_s,R_simulated, simulation_range,new_lmb,R_filtered)
%load('C:\Users\jbojarczuk\Desktop\!testy_do_JLT.mat')
%lambdas = lambdas * 1e9;

%iteracja po kazdym kodzie, ergo - po kazdej lambdzie

Z = zeros(size(dane_po_przetworzeniu,1),length(length_domain));
for i=1:size(dane_po_przetworzeniu,1)
    xcv = varFilter(dane_po_przetworzeniu(i,:),size(dane_referencyjne,2),mean(diff(D_s)),true,true,false);
    xcv_filter = ones(size(xcv));
    xcv_filter(1:D_s(1)-1) = 0;
    xcv = xcv .* xcv_filter;
    Z(i, :) = xcv;
end
%%
Z3 = [dane_po_przetworzeniu-movmean(dane_po_przetworzeniu,D_s(2)-D_s(1),2)];
%% rysowanie z plot3

figure('Renderer', 'opengl', 'units','normalized','Color','w');%,'outerposition',[0 0 1 1]
for i=1:length(length_domain)
    if length_domain(i)>=390 && length_domain(i)<530
        p3 = plot3(repmat(length_domain(i), length(lambdas), 1), lambdas*1e9, Z3(:,i), '*','linewidth', 2, 'color','black');%, 'edgecolor', 'none'
        zlim([0.2 2])
        hold on;
    end
end

%%
%figure;plot(sum(R_simulated{1}.*simulation_range*1e9)/sum(R_simulated{1})-sum(R_received.*lambdas,2)./sum(R_received,2))
R_org = cell2mat(R_simulated');
R_rec = (R_received'-min(R_received'))';
a = sum(R_rec.*lambdas*1e9,2)./sum(R_rec,2)-sum(R_org.*simulation_range*1e9,2)./sum(R_org,2);
figure;
plot(sum(R_rec.*lambdas*1e9,2)./sum(R_rec,2)-sum(R_org.*simulation_range*1e9,2)./sum(R_org,2))
%plot(a([1:12 14:end]))
%%
figure;surf(lambdas*1e9,length_domain(D_s),R_received,'edgecolor','none');
colormap gray

%%
figure('Renderer', 'opengl', 'units','normalized','outerposition',[0 0 1 1],'Color','w');
surf(lambdas*1e9,length_domain,Z','EdgeColor','none');
ylim([390 530])
view([0 90])
shading interp
colormap jet
%%
figure('Renderer', 'opengl', 'units','normalized','outerposition',[0 0 1 1],'Color','w');
surf(lambdas*1e9,length_domain,Z3','EdgeColor','none');
ylim([390 530])
view([0 90])
colormap jet

%%
end

