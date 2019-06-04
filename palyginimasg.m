for alg=1:3
    for i=1:20
        temp=num2str(i,'%03.f')
        switch alg
            case 1
            [ats,laikas(i,alg)]=FAcompare(['dataset\' temp]);
            case 2
            [ats,p,laikas(i,alg)]=myACOcompare(['dataset\' temp]);
            case 3
            [ats,r,laikas(i,alg)]=myMHSAcompare(['dataset\' temp]);
        end
    end
end
mlaik=mean(laikas);
devlaik=std(laikas);
