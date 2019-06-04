for alg=1:3
    for i=1:100
    maxacc=0;
        for j=1:3
            temp=num2str(i,'%03.f')
            switch alg
                case 1
                [ats,laikas(i,alg)]=FAcompare(['dataset\' temp]);
                case 2
                [ats,p,laikas(i,alg)]=myACOcompare(['dataset\' temp]);
                case 3
                [ats,r,laikas(i,alg)]=myMHSAcompare(['dataset\' temp]);
            end
            xe = imread(['dataset\' temp '.png']);
            if(alg==2)
                [Gx,Gy]=imgradientxy(xe);
                Gx=abs(Gx);
                Gy=abs(Gy);
                xe=Gx+Gy;
            end
            TP=0;
            FN=0;
            FP=0;
            for ir=1:120
                for ic=1:120
                    if ats(ir,ic)>0 && xe(ir,ic)>0
                        TP=TP+1;
                        %is the number of pixels that should be included in the segmentation result and are also actually included.
                    elseif ats(ir,ic)>0 && xe(ir,ic)==0
                        FP=FP+1;
                        %is the number of pixels that should be excluded from the segmentation result but actually included
                    end

                end
            end
            ne=sum(sum(ats>0));
            te=sum(sum(xe>0));
            FN=te-TP;
            %is the number of pixels that should be included in the segmentation result but are not
            TN=120*120-FN-TP-FP;
            %True negatives
            accuracy(i,alg)=(TP+TN)/14400; %+% KUR Accuracy?
            if(accuracy(i,alg)>maxacc)
                maxacc=accuracy(i,alg);
                OR(i,alg)=FN/(FN+TP);%+
                UR(i,alg)=FP/(FP+TP);%+
                ER(i,alg)=(FN+FP)/(14400-te);%+
                sensitivity(i,alg)=TP/(TP+FN);%+
                specificity(i,alg)=TN/(TN+FP);%+
                %FPR
                fpr(i,alg)=1-specificity(i,alg);
                npv(i,alg)=TN/(TN+FN);%+
                %RDE Skaiciavimas
                dist=ones(14400)*200;
                %atstumai tarp kiekvieno tasko panasiai kaip ir mhsa
                for ir=1:120
                    for ic=1:120
                        for jr=1:120
                            for jc=1:120
                                if(ats(ir,ic)>0 && xe(jr,jc)>0)
                                    dist(ic+(ir-1)*120,jc+(jr-1)*120)=sqrt((ic-jr)^2+(ir-jr)^2);
                                end
                                %
                                %pazymim taskus kuriuose nera etalono ar segmentuoto vaizdo
                                if(ats(ir,ic)==0||xe(jr,jc)==0)
                                    dist(ic+(ir-1)*120,jc+(jr-1)*120)=200;
                                end
                            end
                        end
                    end
                end
                %minimalus atstumas tarp etalono ir segmentuoto vaizdo tasku
                for temp=1:14400
                    distr(temp)=min(dist(temp,:));
                    distc(temp)=min(dist(:,temp));
                end
                for temp=1:14400
                %atmetam taskus kuriuose nera etalono ar segmentuoto vaizdo
                    if(distc(temp)==200)
                        distc(temp)=0;
                    end
                    if(distr(temp)==200)
                        distr(temp)=0;
                    end
                end
                dr=sum(distr.^2);
                dc=sum(distc.^2);
                RDE(i,alg)=(sqrt(dc/ne)+sqrt(dr/te))/2;
            imwrite(uint8(ats), ['REZ' num2str(alg) 'ITER' num2str(i) '.bmp'], 'bmp');  
            toc
            end
        end
    end
end
mOR=mean(OR);
devOR=std(OR);
mUR=mean(UR);
devUR=std(UR);
mER=mean(ER);
devER=std(ER);
msens=mean(sensitivity);
devSens=std(sensitivity);
mspec=mean(specificity);
devSpec=std(specificity);
macc=mean(accuracy);
devacc=std(accuracy);
mRDE=mean(RDE);
devRDE=std(RDE);
mnpv=mean(npv);
devnpv=std(npv);
mlaik=mean(laikas);
devlaik=std(laikas);