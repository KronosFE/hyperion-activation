#!/usr/bin/env python3
"""HYPERION breeder activation at the QUALIFIED DESIGN-SPEC low-activation alloy (Cr-Ti-V-W, paper 2.48).
Full ENDF/B-VIII.0 chain. 85.04 MW / S_n 3.017e19, wall-loading-anchored flux. 2 full-power-year exposure.
10 CFR 61 §61.55 Class-C classification + long-lived breakout, across nitrogen impurity levels 10/50/100 wppm.
Outputs -> ~/hyperion_designspec_out/
"""
import openmc, openmc.deplete, openmc.data, os, json, csv, numpy as np
np.random.seed(20260726)
OUT=os.path.expanduser("~/hyperion_designspec_out"); os.makedirs(OUT,exist_ok=True); os.chdir(OUT)
openmc.config['cross_sections']=os.environ["OPENMC_CROSS_SECTIONS"]
CHAIN=os.environ["OPENMC_CHAIN"]; openmc.config['chain_file']=CHAIN
try:
    import openmc.deplete.pool as _pool; _pool.USE_MULTIPROCESSING=False
except Exception: pass

# ---- Tier-2 normalization ----
P_FUS_MW=85.04; E_FUS_MEV=17.59; E_N_MEV=14.06; QE=1.602176634e-19
S_n=P_FUS_MW*1e6/(E_FUS_MEV*1e6*QE)      # 3.017e19 n/s
A_fw_cm2=36.06e4
print(f"[norm] S_n={S_n:.4e} n/s  P_n={S_n*E_N_MEV*1e6*QE/1e6:.2f} MW  wall={S_n*E_N_MEV*1e6*QE/1e6/36.06:.3f} MW/m2",flush=True)

# ---- DESIGN alloys (paper 2.48, at%) + achievable impurity spec (wppm) ----
AMU=dict(Cr=52.00,Ti=47.87,V=50.94,W=183.84, Nb=92.91,Mo=95.95,Co=58.93,Ni=58.69,Cu=63.55,Ag=107.87,Zr=91.22,N=14.01)
BASE={"FW":dict(Cr=6.25,Ti=43.75,V=25.0,W=25.0),        # Cr1Ti7V4W4 (first wall / blanket structure)
      "DIV":dict(Cr=12.5,Ti=25.0,V=25.0,W=37.5)}        # Cr2Ti4V4W6 (divertor)
DENS={"FW":7.70,"BLK_front":7.70,"BLK_mid":7.70,"BLK_back":7.70,"DIV":9.5}
IMP_METAL=dict(Nb=5,Mo=50,Co=10,Ni=100,Cu=10,Ag=5,Zr=10)   # wppm, achievable reduced-activation procurement
def mkmat(name, base_key, N_wppm):
    m=openmc.Material(name=name); elems=dict(BASE[base_key])
    Mbar=sum(f*AMU[e] for e,f in elems.items())/sum(elems.values())
    imp=dict(IMP_METAL); imp["N"]=N_wppm
    for e,ppm in imp.items():
        elems[e]=elems.get(e,0.0)+(ppm*1e-6)*Mbar/AMU[e]*100.0
    for e,f in elems.items(): m.add_element(e,f,percent_type='ao')
    m.set_density('g/cm3',DENS[name]); m.depletable=True; return m

def build(N_wppm):
    mfw=mkmat("FW","FW",N_wppm); mb1=mkmat("BLK_front","FW",N_wppm); mb2=mkmat("BLK_mid","FW",N_wppm)
    mb3=mkmat("BLK_back","FW",N_wppm); mdiv=mkmat("DIV","DIV",N_wppm)
    mats=openmc.Materials([mfw,mb1,mb2,mb3,mdiv])
    x0=openmc.XPlane(0.0,boundary_type='vacuum'); xa=openmc.XPlane(2.0); xb=openmc.XPlane(10.0)
    xc=openmc.XPlane(20.0); xd=openmc.XPlane(30.0,boundary_type='vacuum')
    y0=openmc.YPlane(-25,boundary_type='reflective'); y1=openmc.YPlane(25,boundary_type='reflective')
    z0=openmc.ZPlane(-25,boundary_type='reflective'); z1=openmc.ZPlane(25,boundary_type='reflective')
    lat=+y0&-y1&+z0&-z1
    cfw=openmc.Cell(fill=mfw,region=+x0&-xa&lat); c1=openmc.Cell(fill=mb1,region=+xa&-xb&lat)
    c2=openmc.Cell(fill=mb2,region=+xb&-xc&lat); c3=openmc.Cell(fill=mb3,region=+xc&-xd&lat)
    return mats,openmc.Geometry([cfw,c1,c2,c3]),dict(FW=cfw,BLK_front=c1,BLK_mid=c2,BLK_back=c3),\
           dict(FW=mfw,BLK_front=mb1,BLK_mid=mb2,BLK_back=mb3,DIV=mdiv)

# ---- transport once (flux ~ independent of trace impurities) at N=50 baseline ----
mats,geom,cells,matmap=build(50); mats.export_to_xml(); geom.export_to_xml()
src=openmc.IndependentSource(space=openmc.stats.Box((0.01,-25,-25),(0.01,25,25)),
    angle=openmc.stats.Monodirectional((1,0,0)),energy=openmc.stats.Discrete([14.06e6],[1.0]))
st=openmc.Settings(); st.run_mode='fixed source'; st.batches=100; st.particles=int(2e5); st.source=src; st.export_to_xml()
egrid=np.logspace(-2,np.log10(2.0e7),301); ef=openmc.EnergyFilter(egrid); tals=openmc.Tallies()
for nm,c in cells.items():
    t=openmc.Tally(name=f"flux_{nm}"); t.filters=[openmc.CellFilter([c.id]),ef]; t.scores=['flux']; tals.append(t)
tals.export_to_xml(); openmc.run(output=True)
np.save("egrid.npy",egrid); sp=openmc.StatePoint("statepoint.100.h5")
model_area=2500.0; I_src=S_n/A_fw_cm2; S_n_model=I_src*model_area
layers={"FW":(0,2),"BLK_front":(2,10),"BLK_mid":(10,20),"BLK_back":(20,30)}
region_flux={}
for nm in cells:
    fl=sp.get_tally(name=f"flux_{nm}").mean.ravel(); V=(layers[nm][1]-layers[nm][0])*model_area
    phi=fl/V*S_n_model; region_flux[nm]=(phi,float(phi.sum()))
    print(f"[flux] {nm}: {phi.sum():.3e} n/cm2/s",flush=True)
json.dump({"egrid_eV":list(map(float,egrid)),"region_flux_total":{k:v[1] for k,v in region_flux.items()},
           "S_n":S_n,"P_fus_MW":P_FUS_MW},open("flux_actinv.json","w"))

# ---- 10 CFR 61 Class-C (Ci/m3) limits, activated-metal where applicable ----
CLASS_C={'C14':80,'Ni59':220,'Ni63':7000,'Nb94':0.2,'Tc99':3,'I129':0.08,'Sr90':7000,'Cs137':4600}
CLASS_A={'C14':0.8,'Ni59':22,'Ni63':3.5,'Nb94':0.02,'Tc99':0.3,'I129':0.008,'Sr90':0.04,'Cs137':1.0}
LONGLIVED=['C14','Nb94','Ni63','Ni59','Zr93','Mo93','Ag108_m1','Mn53','Tc99']  # >100 yr watch-list
BQ_PER_CI=3.7e10
YR=365.25*86400.0; DAY=86400.0
irr_sub=12; irr=[2*YR/irr_sub]*irr_sub                       # 2 FPY
cool_pts=[DAY,7*DAY,30*DAY,YR,10*YR,100*YR,1000*YR]; cool=[cool_pts[0]]+list(np.diff(cool_pts)); steps=irr+cool
IDX={12:"shutdown",13:"1d",14:"1wk",15:"1mo",16:"1yr",17:"10yr",18:"100yr",19:"1000yr"}
mxf="micro.h5"; phi_fw,_=region_flux["FW"]
micro=openmc.deplete.MicroXS.from_multigroup_flux(energies=egrid,multigroup_flux=phi_fw,chain_file=CHAIN,temperature=294)
micro.to_hdf5(mxf)

def spec_activity_by_nuc(r,mat,idx,mass_g):
    out={}
    for nuc in r[0].index_nuc:
        try: atoms=r.get_atoms(mat,nuc)[1][idx]
        except Exception: continue
        if atoms<=0: continue
        try: hl=openmc.data.half_life(nuc)
        except Exception: hl=None
        if not hl or hl<=0 or not np.isfinite(hl): continue
        out[nuc]=(np.log(2)/hl*atoms)/mass_g   # Bq/g
    return out

results={}
for N_wppm in [10,50,100]:
    _,_,_,mm=build(N_wppm); rows=[]
    for nm,mat in mm.items():
        dens=DENS[nm]; phi,phi_tot=region_flux["FW"] if nm=="DIV" else region_flux[nm]; mat.volume=1.0
        op=openmc.deplete.IndependentOperator(openmc.Materials([mat]),[np.array([1.0])],[micro],
            chain_file=CHAIN,normalization_mode='source-rate',reduce_chain_level=6)
        sr=[phi_tot]*irr_sub+[0.0]*len(cool)
        dpath=f"dep_{nm}_N{N_wppm}.h5"
        openmc.deplete.PredictorIntegrator(op,steps,source_rates=sr,timestep_units='s').integrate(output=False,path=dpath)
        r=openmc.deplete.Results(dpath); mid=list(r[0].index_mat.keys())[0]
        for idx,lab in IDX.items():
            sa=spec_activity_by_nuc(r,mid,idx,dens*1.0)   # Bq/g
            # Ci/m3 concentration = Bq/g * g/cm3 * 1e6 cm3/m3 / Bq_per_Ci
            conc={n:v*dens*1e6/BQ_PER_CI for n,v in sa.items()}
            ci_C=sum(conc.get(n,0)/CLASS_C[n] for n in CLASS_C)
            ci_A=sum(conc.get(n,0)/CLASS_A[n] for n in CLASS_A)
            cls="Class A" if ci_A<1 else ("Class C" if ci_C<1 else "GTCC (>Class C)")
            ll={n:conc.get(n,0) for n in LONGLIVED if conc.get(n,0)>0}
            drivers={n:conc.get(n,0)/CLASS_C[n] for n in CLASS_C if conc.get(n,0)/CLASS_C[n]>1e-4}
            rows.append(dict(N_wppm=N_wppm,region=nm,cooling=lab,classC_index=ci_C,classA_index=ci_A,waste_class=cls,
                             C14_Ci_m3=conc.get('C14',0.0),Nb94_Ci_m3=conc.get('Nb94',0.0),Ni63_Ci_m3=conc.get('Ni63',0.0),
                             Mn53_Ci_m3=conc.get('Mn53',0.0),Zr93_Ci_m3=conc.get('Zr93',0.0),
                             top_classC=";".join(f"{k}:{v:.3f}" for k,v in sorted(drivers.items(),key=lambda kv:-kv[1])[:4])))
        print(f"[N={N_wppm}] {nm}: 100yr ClassC={[x['classC_index'] for x in rows if x['region']==nm and x['cooling']=='100yr'][0]:.3f}",flush=True)
    results[N_wppm]=rows
allrows=[r for N in results for r in results[N]]
with open("classC_by_N.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(allrows[0].keys())); w.writeheader()
    for r in allrows: w.writerow(r)
json.dump(results,open("designspec_results.json","w"),indent=2,default=float)
print("DESIGN-SPEC RUN DONE",flush=True)
