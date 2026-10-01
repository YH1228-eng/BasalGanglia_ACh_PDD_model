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
Amn=0.5;Ams=-0.03;tms=3;
Aln=0.5;Alg=-0.03;tlg=3;
%原连接参数与时间延迟
Ask=1; tsk=2; Agk=1; tgk=2; tmk=2; Amk=0.2; Alk=0.2; tlk=2;%pnn输出参数
Aka=0.1; tka=50; Ake=0.1; tke=50; %输入参数
AAsk=-0.2; tsk=2; AAgk=-0.2; tgk=2; 
%新增加连接参数
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

    % 保存当前值
    Ve_cur = Ve(i);
    Ve1_cur = Ve1(i);
    % 计算延迟索引
    idx_tem = max(1, i - tem/(1000*h));
    Qm_delay = Qm(round(idx_tem));
    % 定义源项（不含Ve和Ve1的部分）
    S = alpha*beta*(Aee*Qe(i) + Aem*Qm_delay + Aei*Qi(i));
    % RK4的4个斜率
    k1_Ve = Ve1_cur;
    k1_Ve1 = S - (alpha+beta)*Ve1_cur - alpha*beta*Ve_cur;
    k2_Ve = Ve1_cur + h/2*k1_Ve1;
    k2_Ve1 = S - (alpha+beta)*(Ve1_cur + h/2*k1_Ve1) - alpha*beta*(Ve_cur + h/2*k1_Ve);
    k3_Ve = Ve1_cur + h/2*k2_Ve1;
    k3_Ve1 = S - (alpha+beta)*(Ve1_cur + h/2*k2_Ve1) - alpha*beta*(Ve_cur + h/2*k2_Ve);
    k4_Ve = Ve1_cur + h*k3_Ve1;
    k4_Ve1 = S - (alpha+beta)*(Ve1_cur + h*k3_Ve1) - alpha*beta*(Ve_cur + h*k3_Ve);
    % RK4更新
    Ve(i+1) = Ve_cur + h/6*(k1_Ve + 2*k2_Ve + 2*k3_Ve + k4_Ve);
    Ve1(i+1) = Ve1_cur + h/6*(k1_Ve1 + 2*k2_Ve1 + 2*k3_Ve1 + k4_Ve1);


    % 保存当前值
    Vi_cur = Vi(i);
    Vi1_cur = Vi1(i);
    % 延迟项
    idx_tim = max(1, i - tim/(1000*h));
    Qm_delay = Qm(round(idx_tim));
    % 定义源项
    S = alpha*beta*(Aie*Qe(i) + Aim*Qm_delay + Aii*Qi(i));
    % RK4的4个斜率
    k1_Vi = Vi1_cur;
    k1_Vi1 = S - (alpha+beta)*Vi1_cur - alpha*beta*Vi_cur;
    k2_Vi = Vi1_cur + h/2*k1_Vi1;
    k2_Vi1 = S - (alpha+beta)*(Vi1_cur + h/2*k1_Vi1) - alpha*beta*(Vi_cur + h/2*k1_Vi);
    k3_Vi = Vi1_cur + h/2*k2_Vi1;
    k3_Vi1 = S - (alpha+beta)*(Vi1_cur + h/2*k2_Vi1) - alpha*beta*(Vi_cur + h/2*k2_Vi);
    k4_Vi = Vi1_cur + h*k3_Vi1;
    k4_Vi1 = S - (alpha+beta)*(Vi1_cur + h*k3_Vi1) - alpha*beta*(Vi_cur + h*k3_Vi);
    % RK4更新
    Vi(i+1) = Vi_cur + h/6*(k1_Vi + 2*k2_Vi + 2*k3_Vi + k4_Vi);
    Vi1(i+1) = Vi1_cur + h/6*(k1_Vi1 + 2*k2_Vi1 + 2*k3_Vi1 + k4_Vi1);

    Qa(i+1)=Qa(i)+h*Qa1(i);
    Qa1(i+1)=Qa1(i)+h*(ra*ra*QQa(i)-2*ra*Qa1(i)-ra*ra*Qa(i));

    % 保存当前值
    Va_cur = Va(i);
    Va1_cur = Va1(i);
    % 延迟项
    idx_tal = max(1, i - tal/(1000*h));
    Ql_delay = Ql(round(idx_tal));
    % 定义源项
    S = alpha*beta*(Aaa*Qa(i) + Aal*Ql_delay + Aab*Qb(i));
    % RK4的4个斜率
    k1_Va = Va1_cur;
    k1_Va1 = S - (alpha+beta)*Va1_cur - alpha*beta*Va_cur;
    k2_Va = Va1_cur + h/2*k1_Va1;
    k2_Va1 = S - (alpha+beta)*(Va1_cur + h/2*k1_Va1) - alpha*beta*(Va_cur + h/2*k1_Va);
    k3_Va = Va1_cur + h/2*k2_Va1;
    k3_Va1 = S - (alpha+beta)*(Va1_cur + h/2*k2_Va1) - alpha*beta*(Va_cur + h/2*k2_Va);
    k4_Va = Va1_cur + h*k3_Va1;
    k4_Va1 = S - (alpha+beta)*(Va1_cur + h*k3_Va1) - alpha*beta*(Va_cur + h*k3_Va);
    % RK4更新
    Va(i+1) = Va_cur + h/6*(k1_Va + 2*k2_Va + 2*k3_Va + k4_Va);
    Va1(i+1) = Va1_cur + h/6*(k1_Va1 + 2*k2_Va1 + 2*k3_Va1 + k4_Va1);

    % 保存当前值
    Vb_cur = Vb(i);
    Vb1_cur = Vb1(i);
    % 延迟项
    idx_tbl = max(1, i - tbl/(1000*h));
    Ql_delay = Ql(round(idx_tbl));
    % 定义源项
    S = alpha*beta*(Aba*Qa(i) + Abl*Ql_delay + Abb*Qb(i));
    % RK4的4个斜率
    k1_Vb = Vb1_cur;
    k1_Vb1 = S - (alpha+beta)*Vb1_cur - alpha*beta*Vb_cur;
    k2_Vb = Vb1_cur + h/2*k1_Vb1;
    k2_Vb1 = S - (alpha+beta)*(Vb1_cur + h/2*k1_Vb1) - alpha*beta*(Vb_cur + h/2*k1_Vb);
    k3_Vb = Vb1_cur + h/2*k2_Vb1;
    k3_Vb1 = S - (alpha+beta)*(Vb1_cur + h/2*k2_Vb1) - alpha*beta*(Vb_cur + h/2*k2_Vb);
    k4_Vb = Vb1_cur + h*k3_Vb1;
    k4_Vb1 = S - (alpha+beta)*(Vb1_cur + h*k3_Vb1) - alpha*beta*(Vb_cur + h*k3_Vb);
    % RK4更新
    Vb(i+1) = Vb_cur + h/6*(k1_Vb + 2*k2_Vb + 2*k3_Vb + k4_Vb);
    Vb1(i+1) = Vb1_cur + h/6*(k1_Vb1 + 2*k2_Vb1 + 2*k3_Vb1 + k4_Vb1);


    % 保存当前值
    Vp_cur = Vp(i);
    Vp1_cur = Vp1(i);
    % 延迟项
    idx_tpe = max(1, i - tpe/(1000*h));
    idx_tpa = max(1, i - tpa/(1000*h));
    Qe_delay = Qe(round(idx_tpe));
    Qa_delay = Qa(round(idx_tpa));
    % 定义源项
    S = alpha*beta*(Ape*Qe_delay + Apa*Qa_delay + App*Qp(i));
    % RK4的4个斜率
    k1_Vp = Vp1_cur;
    k1_Vp1 = S - (alpha+beta)*Vp1_cur - alpha*beta*Vp_cur;
    k2_Vp = Vp1_cur + h/2*k1_Vp1;
    k2_Vp1 = S - (alpha+beta)*(Vp1_cur + h/2*k1_Vp1) - alpha*beta*(Vp_cur + h/2*k1_Vp);
    k3_Vp = Vp1_cur + h/2*k2_Vp1;
    k3_Vp1 = S - (alpha+beta)*(Vp1_cur + h/2*k2_Vp1) - alpha*beta*(Vp_cur + h/2*k2_Vp);
    k4_Vp = Vp1_cur + h*k3_Vp1;
    k4_Vp1 = S - (alpha+beta)*(Vp1_cur + h*k3_Vp1) - alpha*beta*(Vp_cur + h*k3_Vp);
    % RK4更新
    Vp(i+1) = Vp_cur + h/6*(k1_Vp + 2*k2_Vp + 2*k3_Vp + k4_Vp);
    Vp1(i+1) = Vp1_cur + h/6*(k1_Vp1 + 2*k2_Vp1 + 2*k3_Vp1 + k4_Vp1);


    % 保存当前值
    Vc_cur = Vc(i);
    Vc1_cur = Vc1(i);
    % 延迟项
    idx_tce = max(1, i - tce/(1000*h));
    idx_tca = max(1, i - tca/(1000*h));
    Qe_delay = Qe(round(idx_tce));
    Qa_delay = Qa(round(idx_tca));
    % 定义源项
    S = alpha*beta*(Ace*Qe_delay + Aca*Qa_delay + Acc*Qc(i));
    % RK4的4个斜率
    k1_Vc = Vc1_cur;
    k1_Vc1 = S - (alpha+beta)*Vc1_cur - alpha*beta*Vc_cur;
    k2_Vc = Vc1_cur + h/2*k1_Vc1;
    k2_Vc1 = S - (alpha+beta)*(Vc1_cur + h/2*k1_Vc1) - alpha*beta*(Vc_cur + h/2*k1_Vc);
    k3_Vc = Vc1_cur + h/2*k2_Vc1;
    k3_Vc1 = S - (alpha+beta)*(Vc1_cur + h/2*k2_Vc1) - alpha*beta*(Vc_cur + h/2*k2_Vc);
    k4_Vc = Vc1_cur + h*k3_Vc1;
    k4_Vc1 = S - (alpha+beta)*(Vc1_cur + h*k3_Vc1) - alpha*beta*(Vc_cur + h*k3_Vc);
    % RK4更新
    Vc(i+1) = Vc_cur + h/6*(k1_Vc + 2*k2_Vc + 2*k3_Vc + k4_Vc);
    Vc1(i+1) = Vc1_cur + h/6*(k1_Vc1 + 2*k2_Vc1 + 2*k3_Vc1 + k4_Vc1);

    % 保存当前值
    Vs_cur = Vs(i);
    Vs1_cur = Vs1(i);
    % 延迟项
    idx_tsc = max(1, i - tsc/(1000*h));
    idx_tsk = max(1, i - tsk/(1000*h));
    idx_tsp = max(1, i - tsp/(1000*h));
    Qc_delay = Qc(round(idx_tsc));
    Qk_delay = Qk(round(idx_tsk));
    Qp_delay = Qp(round(idx_tsp));
    % 定义源项
    S = alpha*beta*(Asc*Qc_delay + Ask*Qk_delay + AAsk*Qk_delay + Asp*Qp_delay);
    % RK4的4个斜率
    k1_Vs = Vs1_cur;
    k1_Vs1 = S - (alpha+beta)*Vs1_cur - alpha*beta*Vs_cur;
    k2_Vs = Vs1_cur + h/2*k1_Vs1;
    k2_Vs1 = S - (alpha+beta)*(Vs1_cur + h/2*k1_Vs1) - alpha*beta*(Vs_cur + h/2*k1_Vs);
    k3_Vs = Vs1_cur + h/2*k2_Vs1;
    k3_Vs1 = S - (alpha+beta)*(Vs1_cur + h/2*k2_Vs1) - alpha*beta*(Vs_cur + h/2*k2_Vs);
    k4_Vs = Vs1_cur + h*k3_Vs1;
    k4_Vs1 = S - (alpha+beta)*(Vs1_cur + h*k3_Vs1) - alpha*beta*(Vs_cur + h*k3_Vs);
    % RK4更新
    Vs(i+1) = Vs_cur + h/6*(k1_Vs + 2*k2_Vs + 2*k3_Vs + k4_Vs);
    Vs1(i+1) = Vs1_cur + h/6*(k1_Vs1 + 2*k2_Vs1 + 2*k3_Vs1 + k4_Vs1);

   % 保存当前值
Vg_cur = Vg(i);
Vg1_cur = Vg1(i);

% 延迟项
idx_tgc = max(1, i - tgc/(1000*h));
idx_tgk = max(1, i - tgk/(1000*h));
idx_tgp = max(1, i - tgp/(1000*h));
Qc_delay = Qc(round(idx_tgc));
Qk_delay = Qk(round(idx_tgk));
Qp_delay = Qp(round(idx_tgp));

% 定义源项（包含四项：Agc*Qc + Agk*Qk + AAgk*Qk + Agp*Qp）
S = alpha*beta*(Agc*Qc_delay + Agk*Qk_delay + AAgk*Qk_delay + Agp*Qp_delay);

% RK4的4个斜率
k1_Vg = Vg1_cur;
k1_Vg1 = S - (alpha+beta)*Vg1_cur - alpha*beta*Vg_cur;

k2_Vg = Vg1_cur + h/2*k1_Vg1;
k2_Vg1 = S - (alpha+beta)*(Vg1_cur + h/2*k1_Vg1) - alpha*beta*(Vg_cur + h/2*k1_Vg);

k3_Vg = Vg1_cur + h/2*k2_Vg1;
k3_Vg1 = S - (alpha+beta)*(Vg1_cur + h/2*k2_Vg1) - alpha*beta*(Vg_cur + h/2*k2_Vg);

k4_Vg = Vg1_cur + h*k3_Vg1;
k4_Vg1 = S - (alpha+beta)*(Vg1_cur + h*k3_Vg1) - alpha*beta*(Vg_cur + h*k3_Vg);

% RK4更新
Vg(i+1) = Vg_cur + h/6*(k1_Vg + 2*k2_Vg + 2*k3_Vg + k4_Vg);
Vg1(i+1) = Vg1_cur + h/6*(k1_Vg1 + 2*k2_Vg1 + 2*k3_Vg1 + k4_Vg1);


    % 保存当前值
Vm_cur = Vm(i);
Vm1_cur = Vm1(i);

% 延迟项
idx_tms = max(1, i - tms/(1000*h));
idx_tmk = max(1, i - tmk/(1000*h));
Qs_delay = Qs(round(idx_tms));
Qk_delay = Qk(round(idx_tmk));

% 定义源项
S = alpha*beta*(Ams*Qs_delay + Amk*Qk_delay + Amn*Qn(i));

% RK4的4个斜率
k1_Vm = Vm1_cur;
k1_Vm1 = S - (alpha+beta)*Vm1_cur - alpha*beta*Vm_cur;

k2_Vm = Vm1_cur + h/2*k1_Vm1;
k2_Vm1 = S - (alpha+beta)*(Vm1_cur + h/2*k1_Vm1) - alpha*beta*(Vm_cur + h/2*k1_Vm);

k3_Vm = Vm1_cur + h/2*k2_Vm1;
k3_Vm1 = S - (alpha+beta)*(Vm1_cur + h/2*k2_Vm1) - alpha*beta*(Vm_cur + h/2*k2_Vm);

k4_Vm = Vm1_cur + h*k3_Vm1;
k4_Vm1 = S - (alpha+beta)*(Vm1_cur + h*k3_Vm1) - alpha*beta*(Vm_cur + h*k3_Vm);

% RK4更新
Vm(i+1) = Vm_cur + h/6*(k1_Vm + 2*k2_Vm + 2*k3_Vm + k4_Vm);
Vm1(i+1) = Vm1_cur + h/6*(k1_Vm1 + 2*k2_Vm1 + 2*k3_Vm1 + k4_Vm1);


   % 保存当前值
Vl_cur = Vl(i);
Vl1_cur = Vl1(i);

% 延迟项
idx_tlg = max(1, i - tlg/(1000*h));
idx_tlk = max(1, i - tlk/(1000*h));
Qg_delay = Qg(round(idx_tlg));
Qk_delay = Qk(round(idx_tlk));

% 定义源项
S = alpha*beta*(Alg*Qg_delay + Alk*Qk_delay + Aln*Qn(i));

% RK4的4个斜率
k1_Vl = Vl1_cur;
k1_Vl1 = S - (alpha+beta)*Vl1_cur - alpha*beta*Vl_cur;

k2_Vl = Vl1_cur + h/2*k1_Vl1;
k2_Vl1 = S - (alpha+beta)*(Vl1_cur + h/2*k1_Vl1) - alpha*beta*(Vl_cur + h/2*k1_Vl);

k3_Vl = Vl1_cur + h/2*k2_Vl1;
k3_Vl1 = S - (alpha+beta)*(Vl1_cur + h/2*k2_Vl1) - alpha*beta*(Vl_cur + h/2*k2_Vl);

k4_Vl = Vl1_cur + h*k3_Vl1;
k4_Vl1 = S - (alpha+beta)*(Vl1_cur + h*k3_Vl1) - alpha*beta*(Vl_cur + h*k3_Vl);

% RK4更新
Vl(i+1) = Vl_cur + h/6*(k1_Vl + 2*k2_Vl + 2*k3_Vl + k4_Vl);
Vl1(i+1) = Vl1_cur + h/6*(k1_Vl1 + 2*k2_Vl1 + 2*k3_Vl1 + k4_Vl1);

    
    % 保存当前值
Vk_cur = Vk(i);
Vk1_cur = Vk1(i);

% 延迟项
idx_tka = max(1, i - tka/(1000*h));
idx_tke = max(1, i - tke/(1000*h));
Qa_delay = Qa(round(idx_tka));
Qe_delay = Qe(round(idx_tke));

% 定义源项
S = alpha*beta*(Aka*Qa_delay + Ake*Qe_delay);

% RK4的4个斜率
k1_Vk = Vk1_cur;
k1_Vk1 = S - (alpha+beta)*Vk1_cur - alpha*beta*Vk_cur;

k2_Vk = Vk1_cur + h/2*k1_Vk1;
k2_Vk1 = S - (alpha+beta)*(Vk1_cur + h/2*k1_Vk1) - alpha*beta*(Vk_cur + h/2*k1_Vk);

k3_Vk = Vk1_cur + h/2*k2_Vk1;
k3_Vk1 = S - (alpha+beta)*(Vk1_cur + h/2*k2_Vk1) - alpha*beta*(Vk_cur + h/2*k2_Vk);

k4_Vk = Vk1_cur + h*k3_Vk1;
k4_Vk1 = S - (alpha+beta)*(Vk1_cur + h*k3_Vk1) - alpha*beta*(Vk_cur + h*k3_Vk);

% RK4更新
Vk(i+1) = Vk_cur + h/6*(k1_Vk + 2*k2_Vk + 2*k3_Vk + k4_Vk);
Vk1(i+1) = Vk1_cur + h/6*(k1_Vk1 + 2*k2_Vk1 + 2*k3_Vk1 + k4_Vk1);
    end
for i=1:(T-1)/h
    qztime(i)=time(0.14/h+i);
    qzQa(i)=Qa(0.14/h+i);
    qzQe(i)=Qe(0.14/h+i);
    qzVa(i)=Va(0.14/h+i);
    qzVe(i)=Ve(0.14/h+i);
    qzVa1(i)=Va1(0.14/h+i);
    qzVe1(i)=Ve1(0.14/h+i);
end
 
n=4;
[b10,a10]=butter(n/2,[0.1*2/Fs 40*2/Fs],'bandpass');
qqzQa=filter(b10,a10,qzQa);
qqzQe=filter(b10,a10,qzQe);
Nfft=2048;window=1000;noverlap=window/2;
[pxx1,f1]=pwelch(qqzQa,window,noverlap,Nfft,Fs);
[pxx2,f2]=pwelch(qqzQe,window,noverlap,Nfft,Fs);
figure
plot(f1,pxx1,'-hr','linewidth',1)
xlabel('Frequency(Hz)')
ylabel('\phi_a PSD(μV^2)')
% set(gca,'ytick',0:0.5:2.5)
set(gca,'xtick',0:5:40)
set(gca,'xlim',[0 40])
set(gca, 'LooseInset', [0,0,0,0]);
hold on
xline(8,'--k','linewidth',2)
text(4.5, 0.02, 'f=8Hz', 'FontSize', 18);
text(7, 0.02, 'ACh_l_k_e=0.2', 'FontSize', 18);
text(10, 0.05, '|ACh_g_k_i|=0.2', 'FontSize', 18);
figure
plot(f2,pxx2,'-hb','linewidth',1)
xlabel('Frequency(Hz)')
ylabel('\phi_e PSD(μV^2)')
% set(gca,'ytick',0:0.5:2.5)
set(gca,'xtick',0:5:40)
set(gca, 'LooseInset', [0,0,0,0]);
set(gca,'xlim',[0 40])
hold on
xline(8,'--k','linewidth',2)
text(4.5, 0.02, 'f=8Hz', 'FontSize', 18);
text(7, 0.02, 'ACh_l_k_e=0.2', 'FontSize', 18);
text(10, 0.05, '|ACh_g_k_i|=0.2', 'FontSize', 18);
figure
plot(qzVa1,qzVa,'r','linewidth',1.2)
xlabel('dV_a/dt')
ylabel('V_a')
xlim([-2500 2500])
ylim([-25 5])
set(gca, 'LooseInset', [0,0,0,0]);
text(4.5, 3, 'ACh_l_k_e=0.2', 'FontSize', 18);
text(4.5, 2, '|ACh_g_k_i|=0.2', 'FontSize', 18);
figure
plot(qzVe1,qzVe,'b','linewidth',1.2)
xlabel('dV_e/dt')
ylabel('V_e')
xlim([-2500 2500])
ylim([-25 5])
text(4.5, 3, 'ACh_l_k_e=0.2', 'FontSize', 18);
text(4.5, 2, '|ACh_g_k_i|=0.2', 'FontSize', 18);
figure
plot(time,Va,'r','LineWidth',1)
xlabel('time(ms)')
ylabel('V_a(mV)')
xlim([1 10])
ylim([-15 15])
set(gcf,'unit','centimeters','position',[5 5 40 15]);
set(gca, 'LooseInset', [0,0,0,0]);
gtext('ACh_l_k_e=0.2')
gtext('|ACh_g_k_i|=0.2')
ax=axes('Position',[0.1 0.4 0.5 0.2]);
box on;
hold on
plot(time(5000:5500),Va(5000:5500),'r','LineWidth',1)
xlim([5 5.5])
figure
plot(time,Ve,'b','LineWidth',1)
xlabel('time(ms)')
ylabel('V_e(mV)')
xlim([1 10])
ylim([-15 15])
set(gcf,'unit','centimeters','position',[2 2 40 15]);
set(gca, 'LooseInset', [0,0,0,0]);
gtext('ACh_l_k_e=0.2')
gtext('|ACh_g_k_i|=0.2')
ax=axes('Position',[0.1 0.4 0.5 0.2]);
box on;
hold on
plot(time(5000:5500),Ve(5000:5500),'b','LineWidth',1)
xlim([5 5.5])