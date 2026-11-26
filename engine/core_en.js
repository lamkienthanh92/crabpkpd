
// Reference data: 39 CRAB isolates with measured combination MICs (Mero, Col, As = MIC alone; other fields = concentration in each combination)
const DATA = __DATA__;
const DT=0.1, NSTEPS=720, NWIN=240;
const FU_M=0.92, FU_S=0.62, C_SUL=4.0;
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function makeRng(seed){const u=mulberry32(seed);return{u,n(){let a=0;while(a===0)a=u();const b=u();return Math.sqrt(-2*Math.log(a))*Math.cos(2*Math.PI*b);}};}
function lncv(rng,cv){const s=Math.sqrt(Math.log(1+cv*cv));return Math.exp(-0.5*s*s+s*rng.n());}
function lnsig(rng,med,sig){return med*Math.exp(sig*rng.n());}

const MERO={M1:{label:'1 g q8h, 30-min infusion',gd:3,b:[1000,8,0.5]},M2:{label:'1 g q8h, 3-h infusion',gd:3,b:[1000,8,3]},
 M3:{label:'2 g q8h, 3-h infusion',gd:6,b:[2000,8,3]},M4:{label:'2 g q8h, 4-h infusion',gd:6,b:[2000,8,4]},
 M5:{label:'6 g/day, continuous infusion',gd:6,b:[6000,24,24]},M6:{label:'2 g q6h, 3-h infusion',gd:8,b:[2000,6,3]}};
const SUL={S1:{label:'1 g q6h, 30-min infusion',gd:4,b:[1000,6,0.5]},S2:{label:'1 g q6h, 4-h infusion',gd:4,b:[1000,6,4]},
 S3:{label:'3 g q8h, 4-h infusion (9 g/day)',gd:9,b:[3000,8,4]},S4:{label:'2 g q8h, 4-h infusion (6 g/day)',gd:6,b:[2000,8,4]},
 S5:{label:'9 g/day, continuous infusion',gd:9,b:[9000,24,24]},S6:{label:'4 g q8h, 4-h infusion (12 g/day)',gd:12,b:[4000,8,4]}};
const CMS=[0,4.5,6,9,10.8,12];   // 0 = renally adjusted algorithm dose (Nation 2017, Table 3)
const MG_PER_MU=33;
const MERO_IDS=Object.keys(MERO), SUL_IDS=Object.keys(SUL);
function meroReg(id,crcl){const [d,tau,ti]=MERO[id].b;if(crcl>=51)return[d,tau,ti];
  if(crcl>=30){return id==='M5'?[4000,24,24]:[d,12,ti];}
  return id==='M5'?[2000,24,24]:[d/2,12,ti];}
function sulReg(id,crcl){const [d,tau,ti]=SUL[id].b;if(crcl<30){return id==='S5'?[4500,24,24]:[d,12,ti];}return[d,tau,ti];}
function regText(r){const [d,tau,ti]=r;const ds=d>=1000?(d/1000)+' g':d+' mg';
  if(tau===24&&ti===24)return ds+'/day, continuous infusion';
  return ds+' q'+tau+'h, '+(ti<1?'30-min':ti+'-h')+' infusion';}
function cmsLookupMG(c){return c<5?130:c<10?145:c<20?160:c<30?175:c<40?195:c<50?220:c<60?245:c<70?275:c<80?300:c<90?340:360;}  // Nation 2017 Table 3, mg CBA/day for Css 2 mg/L
function cmsEff(mu,crcl){const alg=cmsLookupMG(crcl)/MG_PER_MU;if(mu===0)return alg;return crcl<60?Math.min(mu,alg):mu;}
function colCss(doseMg,cr,z){return doseMg/Math.pow(10,0.0048*cr+1.667)*Math.exp(0.415*z);}   // Nation 2017 regression (Fig 2B) with log-normal scatter (sigma 0.415)

function adjWeight(wt,heightCm,sex){const ibw=(sex==='F'?45.5:50)+2.3*(heightCm/2.54-60);return wt>1.2*ibw?ibw+0.4*(wt-ibw):wt;}   // Crandon 2011: adjusted body weight when actual weight exceeds ideal by more than 20%
function makePatients(n,crcl,crclCv,wt,seed,age,pOld,opts){
  opts=opts||{};const effWt=opts.height?adjWeight(wt,opts.height,opts.sex):wt;
  const rng=makeRng(seed);
  const P={n,crcl:new Float64Array(n),m_k10:new Float64Array(n),m_k12:new Float64Array(n),m_k21:new Float64Array(n),m_V1:new Float64Array(n),
    s_k10:new Float64Array(n),s_k12:new Float64Array(n),s_k21:new Float64Array(n),s_V1:new Float64Array(n),eta:new Float64Array(n),
    rM_elf:new Float64Array(n),rM_csf:new Float64Array(n),rM_tis:new Float64Array(n),rS_elf:new Float64Array(n),old:new Uint8Array(n)};
  const rho=opts.rho||0,sq=Math.sqrt(rho),sr=Math.sqrt(1-rho); const mx=(z,z0)=>rho>0?sq*z0+sr*z:z;
  for(let i=0;i<n;i++){
    const zc=rho>0?rng.n():0,zv=rho>0?rng.n():0;
    const cr=opts.range?(opts.range[0]+(opts.range[1]-opts.range[0])*rng.u()):Math.min(400,Math.max(5,crcl*lncv(rng,crclCv))); P.crcl[i]=cr;
    P.m_V1[i]=0.239*effWt*Math.exp(0.52*mx(rng.n(),zv)); P.m_k10[i]=(0.3922+0.0025*cr)*Math.exp(0.35*mx(rng.n(),zc)); P.m_k12[i]=0.487*Math.exp(0.42*rng.n()); P.m_k21[i]=0.647*Math.exp(0.54*rng.n());
    // Sulbactam: Jaruratanasirikul 2016 (PMID 27671056). CL depends on age group and CrCl; Vc 14.56 L, Vp 9.55 L, Q 11.5 L/h.
    const old=(pOld!=null&&pOld!==undefined)?(rng.u()<pOld):(age>65); P.old[i]=old?1:0;
    const CLs=(old?Math.exp(0.54+0.28*cr/30):Math.exp(1.57+0.46*cr/83))*Math.exp(0.484*mx(rng.n(),zc));
    const V1s=14.56*Math.exp(0.52*mx(rng.n(),zv)),Vp=9.55*Math.exp(0.531*rng.n()),Q=11.5*Math.exp(0.514*rng.n());
    P.s_V1[i]=V1s; P.s_k10[i]=CLs/V1s; P.s_k12[i]=Q/V1s; P.s_k21[i]=Q/Vp;
    P.eta[i]=Math.exp(0.415*mx(rng.n(),zc));
    P.rM_elf[i]=lnsig(rng,0.27,0.6); P.rM_csf[i]=lnsig(rng,0.09,0.65); P.rM_tis[i]=lnsig(rng,0.6,0.4); P.rS_elf[i]=lnsig(rng,0.5,0.4);
  }
  return P;
}
function simulate(n,k10,k12,k21,V1,reg){
  const dose=reg[0],tau=reg[1],tinf=reg[2];
  const ts=Math.round(tau/DT),ti=Math.round(tinf/DT),rate=dose/tinf;
  const out=new Float32Array(n*NWIN);
  for(let i=0;i<n;i++){
    const a=k10[i]+k12[i]+k21[i],disc=Math.sqrt(a*a-4*k10[i]*k21[i]);
    const l1=(-a+disc)/2,l2=(-a-disc)/2;
    const m11=-(k10[i]+k12[i]),m12=k21[i],m21=k12[i],m22=-k21[i];
    const e1=Math.exp(l1*DT),e2=Math.exp(l2*DT),dl=l1-l2;
    const p11=(e1*(m11-l2)-e2*(m11-l1))/dl,p12=(e1*m12-e2*m12)/dl,p21=(e1*m21-e2*m21)/dl,p22=(e1*(m22-l2)-e2*(m22-l1))/dl;
    const det=k10[i]*k21[i],r1=p11-1,r2=p21;
    const g1=(m22*r1-m12*r2)/det,g2=(-m21*r1+m11*r2)/det;
    let x1=0,x2=0;const inv=1/V1[i];const o=i*NWIN;
    for(let s=0;s<NSTEPS;s++){
      const u=((s%ts)<ti)?rate:0;
      const nx1=p11*x1+p12*x2+g1*u,nx2=p21*x1+p22*x2+g2*u;x1=nx1;x2=nx2;
      if(s>=NSTEPS-NWIN)out[o+s-(NSTEPS-NWIN)]=x1*inv;
    }
  }
  return out;
}
function minOf(ser,n){const m=new Float64Array(n);for(let i=0;i<n;i++){let v=1e9;const b=i*NWIN;for(let t=0;t<NWIN;t++){if(ser[b+t]<v)v=ser[b+t];}m[i]=v;}return m;}
function okTime(ser,n,ratio,thr,T){ // ratio: mảng theo bệnh nhân; thr: mảng hoặc số
  const o=new Uint8Array(n);const isArr=typeof thr!=='number';
  for(let i=0;i<n;i++){const th=(isArr?thr[i]:thr)/ratio[i];const b=i*NWIN;let c=0;for(let t=0;t<NWIN;t++){if(ser[b+t]>th)c++;}o[i]=(c/NWIN>=T-1e-9)?1:0;}
  return o;
}
function ratioM(P,site){const n=P.n;if(site==='elf')return P.rM_elf;if(site==='csf')return P.rM_csf;if(site==='tissue')return P.rM_tis;return new Float64Array(n).fill(FU_M);}
function ratioS(P,site,sulCsf){const n=P.n;if(site==='elf')return P.rS_elf;if(site==='csf')return new Float64Array(n).fill(sulCsf);if(site==='tissue')return new Float64Array(n).fill(0.5);return new Float64Array(n).fill(FU_S);}
function sulTarget(site){return site==='elf'?0.475:0.6;}
function cssArr(P,mu,nominal){const a=new Float64Array(P.n),D=cmsEff(mu,nominal)*MG_PER_MU;for(let i=0;i<P.n;i++){const cr=P.crcl[i];a[i]=D/Math.pow(10,0.0048*cr+1.667)*P.eta[i];}return a;}
function okCol(css,site,cc,o){const n=css.length;const out=new Uint8Array(n);
  for(let i=0;i<n;i++){let v;
    if(site==='elf')v=o.rElf*o.fu*css[i];else if(site==='csf')v=0.09*css[i];else if(site==='tissue')v=0.12*o.fu*css[i];else v=o.fu*css[i];
    out[i]=v>=(typeof cc==='number'?cc:cc[i])?1:0;}
  return out;}
function colEff(css,site,o){const n=css.length,a=new Float64Array(n);
  for(let i=0;i<n;i++){a[i]=site==='elf'?o.rElf*o.fu*css[i]:site==='csf'?0.09*css[i]:site==='tissue'?0.12*o.fu*css[i]:o.fu*css[i];}return a;}
function wilson(k,n){const z=1.96,p=k/n,d=1+z*z/n,c=p+z*z/(2*n),h=z*Math.sqrt(p*(1-p)/n+z*z/(4*n*n));return[Math.max(0,(c-h)/d),Math.min(1,(c+h)/d)];}
function mean(u,n){let s=0;for(let i=0;i<n;i++)s+=u[i];return s/n;}
function ccFromMic(micC,mult){return Math.max(0.03125,0.125*micC)*mult;}
// Formulas: which drugs are active and which measured columns give their concentrations in the combination
const FORMS={
 T:{id:'T',name:'Triple combination',short:'Colistin + meropenem + A/S',m:1,s:1,c:1,cmF:'M3',ccF:'C3',caF:'A3'},
 CM:{id:'CM',name:'Colistin + meropenem',short:'Colistin + meropenem',m:1,s:0,c:1,cmF:'MmC',ccF:'CmM'},
 CS:{id:'CS',name:'Colistin + ampicillin/sulbactam',short:'Colistin + A/S',m:0,s:1,c:1,ccF:'CcA',caF:'AcC'},
 MS:{id:'MS',name:'Meropenem + ampicillin/sulbactam',short:'Meropenem + A/S',m:1,s:1,c:0,cmF:'MmA',caF:'AmM'}};
const FORM_IDS=['T','CM','CS','MS'];
const SYN={T:[39,39],CM:[31,39],CS:[5,39],MS:[3,39]};   // isolates with synergy / isolates tested
// Predict combination concentrations from the three monotherapy MICs: kernel on the 39 isolates; meropenem keeps the observed fold reduction, colistin and ampicillin take the neighbour's measured value
function sampleComb(fid,micM,micC,micA,n,rng){
  const F=FORMS[fid],N=DATA.length,lm=Math.log2(micM),lc=Math.log2(micC),la=Math.log2(micA);
  const w=new Float64Array(N);let sw=0,sw2=0;
  for(let j=0;j<N;j++){const d=DATA[j];
    const dm=(lm-Math.log2(d.Mero))/0.35,dc=(lc-Math.log2(d.Col))/1.0,da=(la-Math.log2(d.As))/1.0;
    w[j]=Math.exp(-0.5*(dm*dm+dc*dc+da*da));sw+=w[j];}
  if(sw<=1e-300){for(let j=0;j<N;j++){const d=DATA[j];const dm=(lm-Math.log2(d.Mero)),dc=(lc-Math.log2(d.Col)),da=(la-Math.log2(d.As));w[j]=Math.exp(-0.5*(dm*dm*0.35*0.35+dc*dc+da*da)/1.0);sw+=w[j];}}
  const cum=new Float64Array(N);let acc=0;for(let j=0;j<N;j++){acc+=w[j]/sw;cum[j]=acc;sw2+=(w[j]/sw)*(w[j]/sw);}
  const cm=F.m?new Float64Array(n):null,cc=F.c?new Float64Array(n):null,ca=F.s?new Float64Array(n):null;
  const cscale=Math.max(1,micC/0.25);
  for(let i=0;i<n;i++){const u=rng.u();let j=0;while(j<N-1&&cum[j]<u)j++;const d=DATA[j];
    if(F.m)cm[i]=d[F.cmF]/d.Mero*micM; if(F.c)cc[i]=d[F.ccF]*cscale; if(F.s)ca[i]=d[F.caF];}
  return{cm,cc,ca,neff:1/sw2};
}
function inRange(micM,micC,micA){return micM>=32&&micM<=128&&micC>=0.0625&&micC<=0.25&&micA>=32&&micA<=128;}
function quantile(sorted,q){const k=(sorted.length-1)*q,f=Math.floor(k),c=Math.ceil(k);return sorted[f]+(sorted[c]-sorted[f])*(k-f);}
function bandSeries(ser,n,ratio){const out=[];const buf=new Float32Array(n);
  for(let t=0;t<NWIN;t+=2){for(let i=0;i<n;i++)buf[i]=ser[i*NWIN+t]*ratio[i];const s=Float32Array.from(buf).sort();out.push([t*DT,quantile(s,0.05),quantile(s,0.5),quantile(s,0.95)]);}
  return out;}

function computeAll(inp){
  const f=FORMS[inp.form||'T'];
  const cScale=inp.cScale||1,target=inp.target||0.9,limM=inp.limM!=null?inp.limM:0.2,limC=inp.limC!=null?inp.limC:0.3;
  const n=inp.n,P=makePatients(n,inp.crcl,inp.crclCv/100,inp.weight,20261112,inp.age||65,null,{height:inp.height||165,sex:inp.sex});
  if(inp.sulElfR&&inp.sulElfR!==0.5){const k=inp.sulElfR/0.5;for(let i=0;i<n;i++)P.rS_elf[i]*=k;}
  const rng=makeRng(777),ones=new Uint8Array(n).fill(1);
  let cmA=null,ccA=null,caA=null,neff=null;
  if(inp.mode==='manual'){if(f.m)cmA=new Float64Array(n).fill(inp.cmManual);if(f.c)ccA=new Float64Array(n).fill(inp.ccManual);if(f.s)caA=new Float64Array(n).fill(inp.caManual||8);}
  else{const s=sampleComb(f.id,inp.micM,inp.micC,inp.micA,n,rng);cmA=s.cm;ccA=s.cc;caA=s.ca;neff=s.neff;}
  const ccT=f.c?Float64Array.from(ccA,v=>v*inp.mult):null, csT=f.s?Float64Array.from(caA,v=>v/2):null;
  const site=inp.site,o={fu:inp.fuCol,rElf:inp.colElfR};
  const rM=ratioM(P,site),rS=ratioS(P,site,inp.sulCsf);
  const mSer={},mReg={},mTox={},okM={},sSer={},sReg={},okS={},cssK={},cTox={},okC={};
  for(const id of MERO_IDS){const reg=meroReg(id,inp.crcl);mReg[id]=reg;mSer[id]=simulate(n,P.m_k10,P.m_k12,P.m_k21,P.m_V1,reg);mTox[id]=mean(minOf(mSer[id],n).map(v=>v>=45?1:0),n);okM[id]=f.m?okTime(mSer[id],n,rM,cmA,inp.tgtM):ones;}
  for(const id of SUL_IDS){const reg=sulReg(id,inp.crcl);sReg[id]=reg;sSer[id]=simulate(n,P.s_k10,P.s_k12,P.s_k21,P.s_V1,reg);okS[id]=f.s?okTime(sSer[id],n,rS,csT,sulTarget(site)):ones;}
  for(const mu of CMS){const css=cssArr(P,mu,inp.crcl);if(cScale!==1){for(let i=0;i<css.length;i++)css[i]*=cScale;}cssK[mu]=css;let k=0;for(let i=0;i<n;i++)if(css[i]>=4)k++;cTox[mu]=k/n;okC[mu]=f.c?okCol(css,site,ccT,o):ones;}
  const a=okM[inp.mId],b=okS[inp.sId],c=okC[inp.kMU];let J=0;for(let i=0;i<n;i++)J+=a[i]&b[i]&c[i];J/=n;
  const jci=wilson(Math.round(J*n),n);
  const sel={joint:J,jlo:jci[0],jhi:jci[1],pM:f.m?mean(a,n):null,pS:f.s?mean(b,n):null,pC:f.c?mean(c,n):null,toxM:f.m?mTox[inp.mId]:0,toxC:f.c?cTox[inp.kMU]:0,medCss:(()=>{const s=Float64Array.from(cssK[inp.kMU]).sort();return quantile(s,0.5);})()};
  const combos=[];
  for(const mid of (f.m?MERO_IDS:[null]))for(const sid of (f.s?SUL_IDS:[null]))for(const mu of (f.c?CMS:[null])){
    const x=mid?okM[mid]:ones,y=sid?okS[sid]:ones,z=mu!==null?okC[mu]:ones;let j=0;for(let i=0;i<n;i++)j+=x[i]&y[i]&z[i];j/=n;
    const tm=mid?mTox[mid]:0,tc=mu!==null?cTox[mu]:0;
    combos.push({mid,sid,mu,joint:j,cost:(mid?MERO[mid].gd/6:0)+(sid?SUL[sid].gd/9:0)+(mu!==null?cmsEff(mu,inp.crcl)/9:0),toxM:tm,toxC:tc,safe:(tm<=limM&&tc<limC)});}
  const safeOk=combos.filter(x=>x.safe&&x.joint>=target).sort((p,q)=>p.cost-q.cost||q.joint-p.joint);
  const bestAny=combos.filter(x=>x.safe).sort((p,q)=>q.joint-p.joint||p.cost-q.cost);
  const top5=bestAny.slice(0,5);
  const nS=Math.min(n,600),MM=[16,32,64,128,256],CC=[0.0625,0.125,0.25,0.5,1],AA=[32,64,128,256];
  const ax=f.id==='T'?{r:'M',c:'C',t:'A'}:f.id==='CM'?{r:'M',c:'C',t:null}:f.id==='CS'?{r:'A',c:'C',t:null}:{r:'M',c:'A',t:null};
  const vals={M:MM,C:CC,A:AA},names={M:'Meropenem MIC',C:'Colistin MIC',A:'Ampicillin MIC'};
  const tabs=ax.t?vals[ax.t]:[0];const scan={};
  const serM=mSer[inp.mId],serS=sSer[inp.sId];
  for(const T of tabs){scan[T]=[];for(const R of vals[ax.r]){const row=[];for(const Cc of vals[ax.c]){
    const mic={M:inp.micM,C:inp.micC,A:inp.micA};mic[ax.r]=R;mic[ax.c]=Cc;if(ax.t)mic[ax.t]=T;
    let cmv=null,ccv=null,csv=null;
    if(inp.mode==='manual'){if(f.m)cmv=new Float64Array(nS).fill(inp.cmManual);if(f.c)ccv=new Float64Array(nS).fill(inp.ccManual*inp.mult);if(f.s)csv=new Float64Array(nS).fill((inp.caManual||8)/2);}
    else{const s=sampleComb(f.id,mic.M,mic.C,mic.A,nS,makeRng(1000+Math.round(mic.M*7)+Math.round(mic.C*1000)+mic.A));if(f.m)cmv=s.cm;if(f.c)ccv=Float64Array.from(s.cc,v=>v*inp.mult);if(f.s)csv=Float64Array.from(s.ca,v=>v/2);}
    const om=f.m?okTime(serM,nS,rM,cmv,inp.tgtM):null,os=f.s?okTime(serS,nS,rS,csv,sulTarget(site)):null,oc=f.c?okCol(cssK[inp.kMU].subarray(0,nS),site,ccv,o):null;
    let j=0;for(let i=0;i<nS;i++)j+=(om?om[i]:1)&(os?os[i]:1)&(oc?oc[i]:1);
    const ex=(f.m&&(mic.M<32||mic.M>128))||(f.c&&(mic.C<0.0625||mic.C>0.25))||(f.s&&(mic.A<32||mic.A>128));
    row.push({j:j/nS,ex});}scan[T].push(row);}}
  const map={rows:vals[ax.r],cols:vals[ax.c],tabs:ax.t?tabs:null,rowName:names[ax.r],colName:names[ax.c],tabName:ax.t?names[ax.t]:null,scan,ax};
  const cmp=[],cacheM=new Map(),cacheS=new Map(),cacheC=new Map();
  const gM=t=>{let v=cacheM.get(t);if(!v){v=okTime(mSer[inp.mId],n,rM,t,inp.tgtM);cacheM.set(t,v);}return v;};
  const gS=t=>{let v=cacheS.get(t);if(!v){v=okTime(sSer[inp.sId],n,rS,t,sulTarget(site));cacheS.set(t,v);}return v;};
  const gC=t=>{let v=cacheC.get(t);if(!v){v=okCol(cssK[inp.kMU],site,t,o);cacheC.set(t,v);}return v;};
  for(const fid of FORM_IDS){const F=FORMS[fid];let sum=0,n90=0,sm=0,ss=0,sc=0;
    for(let j=0;j<DATA.length;j++){const d=DATA[j];
      const x=F.m?gM(d[F.cmF]):ones,y=F.s?gS(d[F.caF]/2):ones,z=F.c?gC(d[F.ccF]*inp.mult):ones;let t=0,tm=0,ts=0,tc=0;
      for(let i=0;i<n;i++){t+=x[i]&y[i]&z[i];tm+=x[i];ts+=y[i];tc+=z[i];}
      const p=t/n;sum+=p;if(p>=target)n90++;sm+=tm/n;ss+=ts/n;sc+=tc/n;}
    const N=DATA.length;cmp.push({id:fid,name:F.short,cfr:sum/N,p90:n90/N,pM:F.m?sm/N:null,pS:F.s?ss/N:null,pC:F.c?sc/N:null,syn:SYN[fid]});}
  const chart={mero:f.m?bandSeries(mSer[inp.mId],n,rM):null,sul:f.s?bandSeries(sSer[inp.sId],n,rS):null};
  const med=arr=>{const s=Float64Array.from(arr).sort();return quantile(s,0.5);};
  const cmS=f.m?Float64Array.from(cmA).sort():null;
  return{form:f.id,F:f,cScale,target,limM,limC,colVals:f.c?colEff(cssK[inp.kMU],site,o):null,P,sel,combos,safeOk,bestAny,top5,map,cmp,chart,
    cm:f.m?{p10:quantile(cmS,0.1),p50:quantile(cmS,0.5),p90:quantile(cmS,0.9)}:null,cc:f.c?med(ccT):null,cs:f.s?med(csT):null,neff,
    regM:f.m?mReg[inp.mId]:null,regS:f.s?sReg[inp.sId]:null,effMU:f.c?cmsEff(inp.kMU,inp.crcl):null,
    inRange:inp.mode==='manual'?true:!((f.m&&(inp.micM<32||inp.micM>128))||(f.c&&(inp.micC<0.0625||inp.micC>0.25))||(f.s&&(inp.micA<32||inp.micA>128))),inp};
}
