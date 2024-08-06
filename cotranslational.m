clear all;
close all;
clc;

f1 = fopen('cotranslational.out','w');
fclose(f1);

generationtime = 5400;
max_time = generationtime*10;

mrna_lifetime = 600;
rate_productionpermrna = 0.32; %0.32 for 8000 proteins right before division, 0.16 for 4000 proteins
rate_mrnaproduction = 6/mrna_lifetime;
rate_mrnadecay = 1/mrna_lifetime;
rate_proteindecay = 1/(3600*2);%1/4000;
fusionrate = 1/1000;

%both time50 and tim23 have a halflife of about 2 hours
%Determinants and Regulation of Protein Turnover in Yeast
%by Martin-Perez and Villen


mitosizes = linspace(0.01,0.97,32);
criterialist = 0;%2000;%linspace(500,3000,6);
fusionratelist = [0 logspace(-3.5,-2.5,3)];
%rate_mrnaproduction_list = [6 12 24 48 96]/mrna_lifetime;
rate_productionpermrna_list = [0.1 0.5 1 2 10]*0.32;

%for ci=1:length(rate_productionpermrna_list)
for ci=1:length(fusionratelist)
%for ci=1:length(rate_mrnaproduction_list)
%for ci=1:length(criterialist)
    fusionrate = fusionratelist(ci);
    criteria = 0;%criterialist(ci);
    %rate_mrnaproduction = rate_mrnaproduction_list(ci);
    rate_productionpermrna = rate_productionpermrna_list(ci);
    for si=1:length(mitosizes)
        ci
        si
        
        for ri=1:100

            mito1_size = mitosizes(si);% + rand()*(mitosizes(2) - mitosizes(1));
            mito2_size = 1 - mito1_size;

            mrna_number = 5;
            
            mito1_mrna = 0;
            mito2_mrna = 0;
            for j=1:mrna_number
                if(rand() < mito1_size)
                    mito1_mrna = mito1_mrna + 1;
                else
                    mito2_mrna = mito2_mrna + 1;
                end
            end
            
            initialproteinnumber = 4000; %4000
            mito1_proteins = 0;
            mito2_proteins = 0;
            for pi=1:initialproteinnumber
                if(rand() < mito1_size)
                    mito1_proteins = mito1_proteins + 1;
                else
                    mito2_proteins = mito2_proteins + 1;
                end
            end
            protein_number = mito1_proteins + mito2_proteins;

            count = 1;
            t = 0;
            while(t < max_time)

                if(t > generationtime)
                    mito1_proteins = ceil(mito1_proteins/2);
                    mito2_proteins = ceil(mito2_proteins/2);
                    generationtime = generationtime + 5400;
                end

                total_rate = rate_mrnaproduction + mrna_number*rate_mrnadecay + mrna_number*rate_productionpermrna + protein_number*rate_proteindecay + fusionrate;
                Delta_t = (1/total_rate)*log(1/rand());
                t = t + Delta_t;

                rannum = rand();
                if((rate_mrnaproduction/total_rate) > rannum)
                    %mrna production
                    mrna_number = mrna_number + 1;
                    if(rand() < mito1_size)
                        mito1_mrna = mito1_mrna + 1;
                    else
                        mito2_mrna = mito2_mrna + 1;
                    end
                elseif(((rate_mrnaproduction+mrna_number*rate_mrnadecay)/total_rate) > rannum)
                    %mrna decay
                    mrna_number = mrna_number - 1;
                    if(rand() < (mito1_mrna/(mito1_mrna+mito2_mrna)))
                        mito1_mrna = mito1_mrna - 1;
                    else
                        mito2_mrna = mito2_mrna - 1;
                    end
                elseif(((rate_mrnaproduction+mrna_number*rate_mrnadecay+mrna_number*rate_productionpermrna)/total_rate) > rannum)
                    %protein production
                    %if(rand() < mito1_size)
                    if(rand() < (mito1_mrna/(mito1_mrna+mito2_mrna)))
                        if((mito1_proteins/mito1_size)>criteria)
                            mito1_proteins = mito1_proteins + 1;
                            protein_number = protein_number + 1;
                        end
                    else
                        if((mito2_proteins/mito2_size)>criteria)
                            mito2_proteins = mito2_proteins + 1;
                            protein_number = protein_number + 1;
                        end
                    end
                elseif(((rate_mrnaproduction+mrna_number*rate_mrnadecay+mrna_number*rate_productionpermrna+protein_number*rate_proteindecay)/total_rate) > rannum)
                    %protein decay
                    if(rand() < (mito1_proteins/(mito1_proteins+mito2_proteins)))
                        mito1_proteins = mito1_proteins - 1;
                        protein_number = protein_number - 1;
                    else
                        mito2_proteins = mito2_proteins - 1;
                        protein_number = protein_number - 1;
                    end
                elseif(((rate_mrnaproduction+mrna_number*rate_mrnadecay+mrna_number*rate_productionpermrna+protein_number*rate_proteindecay + fusionrate)/total_rate) > rannum)
                    %fusion
                    if(((mito1_proteins/mito1_size) >= criteria) && ((mito2_proteins/mito2_size) >= criteria))
                        total_proteins = mito1_proteins + mito2_proteins;
                        mito1_proteins = floor(mito1_size*total_proteins + 0.5);
                        mito2_proteins = total_proteins - mito1_proteins;
                    end
                end

        %         t_rec(count) = t;
        %         p1_rec(count) = mito1_proteins;
        %         p2_rec(count) = mito2_proteins;
        %         m_rec(count) = mrna_number;
        %         count = count + 1;
            end

        %     plot(t_rec,p1_rec/mito1_size);
        %     hold on
        %     plot(t_rec,p2_rec/mito2_size);
        % 
        %     figure
        %     plot(t_rec,m_rec);

            p1_record(ri) = mito1_proteins;
            p2_record(ri) = mito2_proteins;

%             plot(mito1_size,(mito1_proteins/mito1_size)/4000,'.k');
%             hold on
%             plot(mito2_size,(mito2_proteins/mito2_size)/4000,'.k');
        end
        %p1_recordkeep(si,:) = p1_record;
        cov(si) = std(p1_record)/mean(p1_record);
        
        cov_rec(ci,si) = cov(si);
    end

    %figure
    plot(mitosizes,cov,'Linewidth',3);
    hold on
end
set(gca,'fontsize',14)
xlabel('Fraction of mitochondrial volume','Fontsize',18);
ylabel('Coefficient of variation of protein number','Fontsize',18);
legend('fusion rate = 0','10^{-3.5}','10^{-3}','10^{-2.5}','Fontsize',14);
%legend('threshold = 500','1000','1500','2000','2500','3000','Fontsize',14);
ylim([0 0.6]);

% figure
% for i=1:32
%     plot(mitosizes,p1_recordkeep(i,:)./mitosizes,'.k');
% end

f1 = fopen('cotranslational.out','a');
for i=1:length(mitosizes)
    fprintf(f1,'%E',mitosizes(i));
    for j=1:length(fusionratelist)
        fprintf(f1,' %E',cov_rec(j,i))
    end
    fprintf(f1,'\n',mitosizes(i));
end
fclose(f1);