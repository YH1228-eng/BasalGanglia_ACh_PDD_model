clc;clear;
T=10;Fs=1000;h=1/Fs;
%步长与时间
alpha=160;beta=640;
Maxe=300;se=14;Maxi=300;si=14;Maxa=300;sa=14;Maxb=300;sb=14;
Maxp=65;sp=19;Maxc=65;sc=19;Maxs=250;ss=10;Maxg=250;sg=10;
Maxm=300;sm=13;Maxl=300;sl=13;v=3.3;
%原最大放电率与斜率
Maxk=200;sk=9;
%新增pnn
re=125;Aee=1.6;tee=0;Aem=0.8;tem=35;Aei=-1.2;tei=0;
Aie=1.6;tie=0;Aim=0.2;tim=35;Aii=-0.2;
ra=125;Aaa=1.6;taa=0;Aal=0.8;tal=35;Aab=-1.2;tab=0;
Aba=1.6;tba=0;Abl=0.2;tbl=35;Abb=-0.2;
Ape=0.1;tpe=2;Apa=0.7;tpa=2;App=-0.3;
Aca=0.1;tca=2;Ace=0.7;tce=2;Acc=-0.3;
Asc=-0.1;tsc=1;Asp=-0.04;tsp=1;
Agp=-0.1;tgp=1;Agc=-0.04;tgc=1;
Amn=0.5;tms=3;%Ams=-0.03;
Aln=0.5;tlg=3;Alg=-0.03;
%原连接参数与时间延迟
Ask=1; tsk=2; Agk=1; tgk=2; tmk=2;  Alk=0.2; tlk=2;%pnn输出参数
Aka=0.1; tka=50; Ake=0.1; tke=50; %输入参数
tsk=2; AAgk=-0.2; tgk=2;
%新增加连接参数
for jj=1:41
    Amk(jj)=0.2-(jj-1)*0.005;
    for j=1:41
        AAsk(j)=0.2-(j-1)*0.005;
        for jjj=1:5
            Ams(jjj)=0.03+(jjj-1)*0.005;
            z=51;
            Ve(1:z)=0;Ve1(1:z)=0;Vi(1:z)=0;Vi1(1:z)=0;
            Va(1:z)=0;Va1(1:z)=0;Vb(1:z)=0;Vb1(1:z)=0;
            Vp(1:z)=0;Vp1(1:z)=0;Vc(1:z)=0;Vc1(1:z)=0;
            Vs(1:z)=0;Vs1(1:z)=0;Vg(1:z)=0;Vg1(1:z)=0;
            Vt(1:z)=0;Vt1(1:z)=0;Vm(1:z)=0;Vm1(1:z)=0;
            Vl(1:z)=0;Vl1(1:z)=0;Qa(1:z)=0;Qa1(1:z)=0;
            Qe(1:z)=0;Qe1(1:z)=0;time(1:z)=0;ti(1:z)=0;
            Vk(1:z)=0;Vk1(1:z)=0;
            %初始值
            for i=1:T/h
                Qn(i)=10+sqrt(0)*randn(1);
            end

            for i=z:T/h
                time(i+1)=h*i;
                ti(i)=h*(i-1);
                QQe(i)=Maxe/(1+exp(((se-Ve(i))/v)));
                Qi(i)=Maxi/(1+exp(((si-Vi(i))/v)));
                QQa(i)=Maxa/(1+exp(((sa-Va(i))/v)));
                Qb(i)=Maxb/(1+exp(((sb-Vb(i))/v)));
                Qp(i)=Maxp/(1+exp(((sp-Vp(i))/v)));
                Qc(i)=Maxc/(1+exp(((sc-Vc(i))/v)));
                Qs(i)=Maxs/(1+exp(((ss-Vs(i))/v)));
                Qg(i)=Maxg/(1+exp(((sg-Vg(i))/v)));
                Qm(i)=Maxm/(1+exp(((sm-Vm(i))/v)));
                Ql(i)=Maxl/(1+exp(((sl-Vl(i))/v)));
                Qk(i)=Maxk/(1+exp(((sk-Vk(i))/v)));
                %最大放电率
                Qe(i+1)=Qe(i)+h*Qe1(i);
                Qe1(i+1)=Qe1(i)+h*(re*re*QQe(i)-2*re*Qe1(i)-re*re*Qe(i));

                Ve(i+1)=Ve(i)+h*Ve1(i);
                Ve1(i+1)=Ve1(i)+h*(alpha*beta*(Aee*Qe(i)+Aem*Qm(i-tem/(1000*h))+Aei*Qi(i))-(alpha+beta)*Ve1(i)-alpha*beta*Ve(i));

                Vi(i+1)=Vi(i)+h*Vi1(i);
                Vi1(i+1)=Vi1(i)+h*(alpha*beta*(Aie*Qe(i)+Aim*Qm(i-tim/(1000*h))+Aii*Qi(i))-(alpha+beta)*Vi1(i)-alpha*beta*Vi(i));

                Qa(i+1)=Qa(i)+h*Qa1(i);
                Qa1(i+1)=Qa1(i)+h*(ra*ra*QQa(i)-2*ra*Qa1(i)-ra*ra*Qa(i));

                Va(i+1)=Va(i)+h*Va1(i);
                Va1(i+1)=Va1(i)+h*(alpha*beta*(Aaa*Qa(i)+Aal*Ql(i-tal/(1000*h))+Aab*Qb(i))-(alpha+beta)*Va1(i)-alpha*beta*Va(i));

                Vb(i+1)=Vb(i)+h*Vb1(i);
                Vb1(i+1)=Vb1(i)+h*(alpha*beta*(Aba*Qa(i)+Abl*Ql(i-tbl/(1000*h))+Abb*Qb(i))-(alpha+beta)*Vb1(i)-alpha*beta*Vb(i));

                Vp(i+1)=Vp(i)+h*Vp1(i);
                Vp1(i+1)=Vp1(i)+h*(alpha*beta*(Ape*Qe(i-tpe/(1000*h))+Apa*Qa(i-tpa/(1000*h))+App*Qp(i))-(alpha+beta)*Vp1(i)-alpha*beta*Vp(i));

                Vc(i+1)=Vc(i)+h*Vc1(i);
                Vc1(i+1)=Vc1(i)+h*(alpha*beta*(Ace*Qe(i-tce/(1000*h))+Aca*Qa(i-tca/(1000*h))+Acc*Qc(i))-(alpha+beta)*Vc1(i)-alpha*beta*Vc(i));

                Vs(i+1)=Vs(i)+h*Vs1(i);
                Vs1(i+1)=Vs1(i)+h*(alpha*beta*(Asc*Qc(i-tsc/(1000*h))+Ask*Qk(i-tsk/(1000*h))-AAsk(j)*Qk(i-tsk/(1000*h))+Asp*Qp(i-tsp/(1000*h)))-(alpha+beta)*Vs1(i)-alpha*beta*Vs(i));

                Vg(i+1)=Vg(i)+h*Vg1(i);
                Vg1(i+1)=Vg1(i)+h*(alpha*beta*(Agc*Qc(i-tgc/(1000*h))+Agk*Qk(i-tgk/(1000*h))+AAgk*Qk(i-tgk/(1000*h))+Agp*Qp(i-tgp/(1000*h)))-(alpha+beta)*Vg1(i)-alpha*beta*Vg(i));

                Vm(i+1)=Vm(i)+h*Vm1(i);
                Vm1(i+1)=Vm1(i)+h*(alpha*beta*(-Ams(jjj)*Qs(i-tms/(1000*h))+Amk(jj)*Qk(i-tmk/(1000*h))+Amn*Qn(i))-(alpha+beta)*Vm1(i)-alpha*beta*Vm(i));

                Vl(i+1)=Vl(i)+h*Vl1(i);
                Vl1(i+1)=Vl1(i)+h*(alpha*beta*(Alg*Qg(i-tlg/(1000*h))+Alk*Qk(i-tlk/(1000*h))+Aln*Qn(i))-(alpha+beta)*Vl1(i)-alpha*beta*Vl(i));

                Vk(i+1)=Vk(i)+h*Vk1(i);
                Vk1(i+1)=Vk1(i)+h*(alpha*beta*(Aka*Qa(i-tka/(1000*h))+Ake*Qe(i-tke/(1000*h)))-(alpha+beta)*Vk1(i)-alpha*beta*Vk(i));
            end
            for i=1:(T-1)/h
                qztime(i)=time(1/h+i);
                qzQa(i)=Qa(1/h+i);
                qzQe(i)=Qe(1/h+i);
                qzVa(i)=Va(1/h+i);
                qzVe(i)=Ve(1/h+i);
            end
            %极值
            Vem(jj,j,jjj)=max(qzVe);
            Vei(jj,j,jjj)=min(qzVe);
            Vam(jj,j,jjj)=max(qzVa);
            Vai(jj,j,jjj)=min(qzVa);
            Amp_a(jj,j,jjj)=Vam(jj,j)-Vai(jj,j);
            Amp_e(jj,j,jjj)=Vem(jj,j)-Vei(jj,j);
            %功率比
            n=4;
            [b10,a10]=butter(n/2,[0.1*2/Fs 40*2/Fs],'bandpass');
            qqzQa=filter(b10,a10,qzQa);
            qqzQe=filter(b10,a10,qzQe);
            Nfft=1024;window=1000;noverlap=window/2;
            [pxx1,f1]=pwelch(qqzQa,window,noverlap,Nfft,Fs);
            [pxx2,f2]=pwelch(qqzQe,window,noverlap,Nfft,Fs);
            g1=0;g2=0;g3=0;g4=0;
            s1=size(pxx1);
            s2=size(pxx2);
            for k=1:s1(1)
                if f1(k)>=1
                    k1=k;
                    break
                end
            end
            for k=k1:s1(1)
                if f1(k)>7.5
                    k2=k;
                    break
                end
            end
            for k=k1:k2-1
                g1=g1+pxx1(k);
            end

            for k=1:s1(1)
                if f1(k)>=8
                    k3=k;
                    break
                end
            end
            for k=k3:s1(1)
                if f1(k)>35.5
                    k4=k;
                    break
                end
            end
            for k=k3:k4-1
                g2=g2+pxx1(k);
            end
            for k=1:s2(1)
                if f2(k)>=1
                    k5=k;
                    break
                end
            end
            for k=k5:s2(1)
                if f2(k)>7.5
                    k6=k;
                    break
                end
            end
            for k=k5:k6-1
                g3=g3+pxx2(k);
            end

            for k=1:s2(1)
                if f2(k)>=8
                    k7=k;
                    break
                end
            end
            for k=k7:s2(1)
                if f2(k)>35.5
                    k8=k;
                    break
                end
            end
            for k=k7:k8-1
                g4=g4+pxx2(k);
            end
            g5(jj,j,jjj)=g2/g1;
            g6(jj,j,jjj)=g4/g3;
        end
    end
end

[Amk,AAsk,Ams]=meshgrid(Amk,AAsk,Ams);
xslice = [];   
yslice = [];
zslice = [0.03, 0.035, 0.04, 0.045, 0.05];
figure
slice(Amk,AAsk,Ams,g5,xslice,yslice,zslice,'nearest')
xlabel('|ACh_s_k_i|(mVs)')
ylabel('ACh_m_k_e(mVs)')
zlabel('|A_m_s|(mVs)')
title('Power ratio of \phi_a')
grid off
figure
slice(Amk,AAsk,Ams,g6,xslice,yslice,zslice,'nearest')
xlabel('|ACh_s_k_i|(mVs)')
ylabel('ACh_m_k_e(mVs)')
zlabel('|A_m_s|(mVs)')
title('Power ratio of \phi_e')
grid off