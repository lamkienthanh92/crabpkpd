// Reproducible comparison of model output with published probabilities of target attainment (Supplementary Table S2).
// Usage: node s2_verification.js > s2_results.json   (run in the folder holding core_en.js and data_js.txt)
const fs=require('fs');
const core=fs.readFileSync('./core_en.js','utf8').replace('__DATA__',fs.readFileSync('./data_js.txt','utf8'));
const C=new Function(core+';return {makePatients,simulate,okTime,mean,FU_M,adjWeight,makeRng,cmsEff,cssArr};')();
const N=10000, ones=v=>new Float64Array(N).fill(v);
function pats(lo,hi,seed,ag,wt){ // wt: fixed weight (kg) or null for the study distribution (71 ± 22 kg)
  const P=C.makePatients(N,(lo+hi)/2,0,1,seed,null,ag,{range:[lo,hi]});
  const rr=C.makeRng(seed+1);
  for(let i=0;i<N;i++){let w,h,male;
    if(wt){P.m_V1[i]*=wt;} else {w=71+22*rr.n();w=Math.min(150,Math.max(35,w));male=rr.u()<0.65;h=165+8*rr.n();P.m_V1[i]*=C.adjWeight(w,Math.min(200,Math.max(140,h)),male?'M':'F');}}
  return P;}
const mero=(lo,hi,reg,mic,T,wt,seed)=>{const P=pats(lo,hi,seed,0,wt);const s=C.simulate(N,P.m_k10,P.m_k12,P.m_k21,P.m_V1,reg);return +(100*C.mean(C.okTime(s,N,ones(C.FU_M),mic,T),N)).toFixed(1);};
const sulPlasma=(lo,hi,reg,mic,T,fu,seed)=>{let a=0;for(const ag of [0,1]){const P=pats(lo,hi,seed+ag,ag,null);const s=C.simulate(N,P.s_k10,P.s_k12,P.s_k21,P.s_V1,reg);a+=C.mean(C.okTime(s,N,ones(fu),mic,T),N);}return +(50*a).toFixed(1);};
const sulElf=(lo,hi,reg,mic,T,ag,seed)=>{const P=pats(lo,hi,seed,ag,null);const s=C.simulate(N,P.s_k10,P.s_k12,P.s_k21,P.s_V1,reg);return +(100*C.mean(C.okTime(s,N,P.rS_elf,mic,T),N)).toFixed(1);};
const sulElfMix=(lo,hi,reg,mic,T,seed)=>+((sulElf(lo,hi,reg,mic,T,0,seed)+sulElf(lo,hi,reg,mic,T,1,seed+7))/2).toFixed(1);
function colShare(lo,hi,doseFn,thr,seed){const P=pats(lo,hi,seed,0,null);let k=0;for(let i=0;i<N;i++){const cr=P.crcl[i];const css=doseFn(cr)/Math.pow(10,0.0048*cr+1.667)*P.eta[i];if(css>=thr)k++;}return +(100*k/N).toFixed(1);}
const alg=cr=>(cr<5?130:cr<10?145:cr<20?160:cr<30?175:cr<40?195:cr<50?220:cr<60?245:cr<70?275:cr<80?300:cr<90?340:360);
const R={};
// Crandon 2011 (40% fT>MIC; CrCl ≥50 modelled as 50–200, alternative 50–120)
const crd=[['2 g q8h, 3 h, CrCl ≥50, MIC 8',[2000,8,3],8,[50,200],[50,120]],['2 g q8h, 3 h, CrCl ≥50, MIC 16',[2000,8,3],16,[50,200],[50,120]],['2 g q8h, 0.5 h, CrCl ≥50, MIC 8',[2000,8,0.5],8,[50,200],[50,120]],
 ['1 g q8h, 3 h, CrCl 30–49, MIC 8',[1000,8,3],8,[30,49],null],['1 g q8h, 0.5 h, CrCl 30–49, MIC 8',[1000,8,0.5],8,[30,49],null],['1 g q12h, 3 h, CrCl 10–29, MIC 8',[1000,12,3],8,[10,29],null],['1 g q12h, 0.5 h, CrCl 10–29, MIC 8',[1000,12,0.5],8,[10,29],null],['1 g q12h, 0.5 h, CrCl 10–29, MIC 4',[1000,12,0.5],4,[10,29],null]];
R.crandon=crd.map(([l,reg,mic,r1,r2],i)=>[l,mero(r1[0],r1[1],reg,mic,0.4,null,100+i),r2?mero(r2[0],r2[1],reg,mic,0.4,null,100+i):null]);
// Luo 2025 (100% fT>MIC; fixed CrCl; weight not stated: 58 and 70 kg)
R.luo=[];for(const wt of [58,70])for(const cr of [140,160,180,200])for(const mic of [1,2,4])R.luo.push([wt,cr,mic,mero(cr,cr,[2000,8,3],mic,1.0,wt,300+cr)]);
// Jaruratanasirikul 2016 (60% T>MIC, plasma; CrCl 50–120, both age groups equally)
R.sulPlasma=[];for(const [l,reg,mic] of [['1 g q6h, 4 h, MIC 4',[1000,6,4],4],['1 g q6h, 4 h, MIC 8',[1000,6,4],8],['1 g q6h, 4 h, MIC 16',[1000,6,4],16],['3 g q8h, 4 h, MIC 32',[3000,8,4],32]])R.sulPlasma.push([l,sulPlasma(50,120,reg,mic,0.6,0.62,500),sulPlasma(50,120,reg,mic,0.6,1.0,500)]);
// Abouelhassan 2024 (ELF; published PK from healthy volunteers, clearance 13.7 L/h)
R.sulElf=[];for(const [l,reg,mic,T] of [['1 g q6h, 0.5 h, MIC 4, target 16.62%',[1000,6,0.5],4,0.1662],['3 g q8h, 4 h, MIC 16, target 47.51%',[3000,8,4],16,0.4751]])R.sulElf.push([l,sulElf(130,180,reg,mic,T,0,700),sulElfMix(50,120,reg,mic,T,710)]);
// Nation 2017 (Css ≥2 mg/L)
R.col=[['Algorithm dose, CrCl 50 to <80',colShare(50,79.99,alg,2,900)],['360 mg CBA/day, CrCl 80–236',colShare(80,236,()=>360,2,910)]];
console.log(JSON.stringify(R));
