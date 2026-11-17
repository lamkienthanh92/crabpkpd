const fs=require('fs');
const core=fs.readFileSync('./core_en.js','utf8').replace('__DATA__',fs.readFileSync('./data_js.txt','utf8'));
const C=new Function(core+';return {makePatients,simulate,meroReg,sulReg,okTime,ratioM,ratioS,sulTarget,cssArr,okCol,cmsEff,adjWeight,MERO,SUL,CMS,MERO_IDS,SUL_IDS,mean,minOf,makeRng,regText,NWIN};')();
const ISO=JSON.parse(fs.readFileSync(process.env.ISO||'iso39.json','utf8')); const NI=ISO.length;
const N=1500;
const STRATA=[['A: CrCl 130–180',130,180,155],['B: CrCl 90–129',90,129,110],['C: CrCl 51–89',51,89,70],['D: CrCl 30–50',30,50,40],['E: CrCl 10–29',10,29,20]];
const SITES=['plasma','elf','csf','tissue'];
const SCEN=(process.env.SCEN||'S0:0.49:1,S3:0.085:4').split(',').map(x=>{const q=x.split(':');return [q[0],+q[1],+q[2]];});
const MERO_T=+(process.env.MERO_T||0.4), ELF_K=+(process.env.ELF_K||1), CSF_K=+(process.env.CSF_K||1), IND_CMS=+(process.env.IND_CMS||0);
// formulas: thresholds per isolate (meropenem ug/mL, sulbactam = ampicillin/2, colistin ug/mL)
const FORM=[
 {id:'MONO_M',name:'Meropenem alone',m:i=>i.Mero,s:null,c:null},
 {id:'MONO_C',name:'Colistin alone',m:null,s:null,c:i=>i.Col},
 {id:'MONO_S',name:'Ampicillin/sulbactam alone',m:null,s:i=>i.As/2,c:null},
 {id:'CM',name:'Colistin + meropenem',m:i=>i.MmC,s:null,c:i=>i.CmM},
 {id:'CS',name:'Colistin + ampicillin/sulbactam',m:null,s:i=>i.AcC/2,c:i=>i.CcA},
 {id:'MS',name:'Meropenem + ampicillin/sulbactam',m:i=>i.MmA,s:i=>i.AmM/2,c:null},
 {id:'T',name:'Colistin + meropenem + ampicillin/sulbactam',m:i=>i.M3,s:i=>i.A3/2,c:i=>i.C3}];
const ones=new Uint8Array(N).fill(1);
const out=[]; const PP={}; const PI={}; const t0=Date.now();
for(let si=0;si<STRATA.length;si++){
  const [sname,lo,hi,nom]=STRATA[si]; const seed=777+si*13;
  const rr=C.makeRng(seed+1); const effW=new Float64Array(N);
  for(let i=0;i<N;i++){let w=71+22*rr.n();w=Math.min(150,Math.max(35,w));const male=rr.u()<0.65;const h=165+8*rr.n();effW[i]=C.adjWeight(w,Math.min(200,Math.max(140,h)),male?'M':'F');}
  const Pb=C.makePatients(N,nom,0,1,seed,null,0,{range:[lo,hi]}); for(let i=0;i<N;i++){Pb.m_V1[i]*=effW[i];Pb.rM_elf[i]*=ELF_K;Pb.rM_csf[i]*=CSF_K;}
  const serM={},toxM={}; for(const id of C.MERO_IDS){const reg=C.meroReg(id,nom);serM[id]=C.simulate(N,Pb.m_k10,Pb.m_k12,Pb.m_k21,Pb.m_V1,reg);toxM[id]=C.mean(C.minOf(serM[id],N).map(v=>v>=45?1:0),N);}
  const cacheM=new Map(); const okMf=(id,site,thr)=>{const k=id+'|'+site+'|'+thr;let v=cacheM.get(k);if(!v){v=C.okTime(serM[id],N,C.ratioM(Pb,site),thr,MERO_T);cacheM.set(k,v);}return v;};
  for(const ag of [0,1]){
    const P=C.makePatients(N,nom,0,1,seed,null,ag,{range:[lo,hi]}); for(let i=0;i<N;i++){P.m_V1[i]*=effW[i];P.rM_elf[i]*=ELF_K;P.rM_csf[i]*=CSF_K;}
    const serS={}; for(const id of C.SUL_IDS){serS[id]=C.simulate(N,P.s_k10,P.s_k12,P.s_k21,P.s_V1,C.sulReg(id,nom));}
    const cacheS=new Map(); const okSf=(id,site,thr)=>{const k=id+'|'+site+'|'+thr;let v=cacheS.get(k);if(!v){v=C.okTime(serS[id],N,C.ratioS(P,site,0.5),thr,C.sulTarget(site));cacheS.set(k,v);}return v;};
    const cssA={},toxC={}; for(const mu of C.CMS){let css=C.cssArr(P,mu,nom); if(IND_CMS){css=new Float64Array(N);for(let i=0;i<N;i++){css[i]=C.cmsEff(mu,P.crcl[i])*33/Math.pow(10,0.0048*P.crcl[i]+1.667)*P.eta[i];}} cssA[mu]=css;let k=0;for(let i=0;i<N;i++)if(css[i]>=4)k++;toxC[mu]=k/N;}
    const cacheC=new Map(); const okCf=(mu,site,cc,fu)=>{const k=mu+'|'+site+'|'+cc+'|'+fu;let v=cacheC.get(k);if(!v){v=C.okCol(cssA[mu],site,cc,{fu:fu,rElf:0.26});cacheC.set(k,v);}return v;};
    for(const F of FORM)for(const site of SITES)for(const [sc,fu,mult] of SCEN){
      const ML=F.m?C.MERO_IDS:[null], SL=F.s?C.SUL_IDS:[null], KL=F.c?C.CMS:[null];
      for(const m of ML)for(const s of SL)for(const k of KL){
        let sum=0,n90=0;
        const pkM={M1:1,M2:2,M4:3}[m],pkS={S1:1,S2:2,S3:3}[s];
        const pk=(m&&s)?(pkM&&pkM===pkS?pkM:0):(m?pkM||0:(s?pkS||0:0));
        const want=(site==='plasma'&&sc==='S0'&&pk>0&&(k===null||k===0));
        const pat=want?new Float64Array(N):null; const piso=want?new Float64Array(NI):null;
        for(let j=0;j<NI;j++){const I=ISO[j];
          const a=F.m?okMf(m,site,F.m(I)):ones, b=F.s?okSf(s,site,F.s(I)):ones, c=F.c?okCf(k,site,F.c(I)*mult,fu):ones;
          let t=0;for(let i=0;i<N;i++){const v=a[i]&b[i]&c[i];t+=v;if(want)pat[i]+=v/NI;} const p=t/N; sum+=p; if(want)piso[j]=p; if(p>=0.9)n90++;}
        if(want){PP[F.id+'|'+sname+'|'+(ag?'>65':'<=65')+'|P'+pk]=Array.from(pat,v=>+v.toFixed(5));PI[F.id+'|'+sname+'|'+(ag?'>65':'<=65')+'|P'+pk]=Array.from(piso,v=>+v.toFixed(5));}
        const tm=m?toxM[m]:0, tc=k!==null?toxC[k]:0;
        out.push([F.id,sname,ag?'>65':'<=65',site,sc,m,s,k,+(sum/NI).toFixed(4),+(n90/NI).toFixed(4),+tm.toFixed(4),+tc.toFixed(4),(tm<=0.2&&tc<0.3)?1:0,
          (m?C.MERO[m].gd/6:0)+(s?C.SUL[s].gd/9:0)+(k!==null?C.cmsEff(k,nom)/9:0)]);
      }
    }
  }
  console.log(sname,'done',Math.round((Date.now()-t0)/1000),'s');
}
fs.writeFileSync(process.env.OUT||'results_formulas.json',JSON.stringify({pp:PP,pi:PI,cols:['formula','stratum','age','site','scen','m','s','k','cfr','p90','toxM','toxC','safe','cost'],rows:out,forms:FORM.map(f=>[f.id,f.name]),N,NI}));
console.log('rows',out.length);
