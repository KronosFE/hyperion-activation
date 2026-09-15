#!/usr/bin/env python3
"""Shutdown-dose back-end for the HYPERION breeder activation package.
Reads the EXISTING per-region depletion inventories (results/dep_*_N50.h5), derives the decay-photon
(gamma) source per region per cooling time, and estimates a CONTACT / surface dose rate via a stated
semi-infinite self-dose model. Adds nothing to and changes nothing in the activation/waste result.

Contact-dose model (stated first-order estimate, NOT a 3-D transport map):
  A component is treated as a semi-infinite homogeneous slab with a uniform volumetric decay-photon
  source S_v(E) [photons/s/cm3]. The surface scalar flux is phi_s(E) = S_v(E) / (2 mu(E)), where mu(E)
  is the material's total photon attenuation coefficient [1/cm] (self-shielding). The contact ambient
  dose rate is H = sum_E phi_s(E) * h(E), with h(E) the ICRP-74 photon fluence-to-ambient-dose-equivalent
  coefficient (AP) [pSv*cm2]. This is the standard R2S contact/self-dose approximation.
"""
import openmc, openmc.deplete, openmc.data, os, json, csv, numpy as np
from collections import defaultdict
PKG=os.environ.get("PKG_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES=os.path.join(PKG,"results")
openmc.config['cross_sections']=os.environ["OPENMC_CROSS_SECTIONS"]
openmc.config['chain_file']=os.environ["OPENMC_CHAIN"]
DENS={"FW":7.70,"BLK_front":7.70,"BLK_mid":7.70,"BLK_back":7.70,"DIV":9.5}
IDX={12:"shutdown",13:"1d",14:"1wk",15:"1mo",16:"1yr",17:"10yr",18:"100yr",19:"1000yr"}

# ICRP-74 photon flux-to-ambient-dose (AP), pSv*cm2
de,dc=openmc.data.dose_coefficients('photon',geometry='AP')
de=np.asarray(de); dc=np.asarray(dc)

# photon attenuation mu(E) [1/cm] from element photon data (total: coherent+incoherent+photoelectric+pair)
LIB=openmc.data.DataLibrary.from_xml(os.environ["OPENMC_CROSS_SECTIONS"]); _pcache={}
ATT_MT=[502,504,515,516,517,522]
def _photon(el):
    if el not in _pcache:
        e=LIB.get_by_material(el,data_type='photon'); _pcache[el]=openmc.data.IncidentPhoton.from_hdf5(e['path'])
    return _pcache[el]
def mu_of(mat,E):   # E in eV, returns 1/cm
    eld=defaultdict(float)
    for nuc,ad in mat.get_nuclide_atom_densities().items():   # atoms/b-cm
        z=openmc.data.zam(nuc)[0]; eld[openmc.data.ATOMIC_SYMBOL[z]]+=ad
    mu=0.0
    for el,ad in eld.items():
        try: ph=_photon(el)
        except Exception: continue
        sig=0.0
        for mt in ATT_MT:
            if mt in ph.reactions:
                try: sig+=float(ph[mt].xs['0K'](E)) if hasattr(ph[mt].xs,'__getitem__') else float(ph[mt].xs(E))
                except Exception: pass
        mu+=ad*sig    # atoms/b-cm * b = 1/cm
    return max(mu,1e-6)

rows=[]; block={}
for nm in DENS:
    r=openmc.deplete.Results(os.path.join(RES,f"dep_{nm}_N50.h5")); block[nm]={}
    mat_id=list(r[0].index_mat.keys())[0]; nucs=list(r[0].index_nuc)
    atoms_ts={nuc:r.get_atoms(mat_id,nuc)[1] for nuc in nucs}   # atoms in the 1 cm3 basis, vs time
    for idx,lab in IDX.items():
        at={nuc:atoms_ts[nuc][idx] for nuc in nucs if atoms_ts[nuc][idx]>0}
        tot=sum(at.values())
        mt=openmc.Material()
        for nuc,a in at.items():
            try: mt.add_nuclide(nuc,a,'ao')
            except Exception: pass
        mt.set_density('atom/b-cm', tot*1e-24); mt.volume=1.0     # atoms/cm3 -> atoms/b-cm
        ph=mt.get_decay_photon_energy()      # photons/s per cm3 (specific basis)
        if ph is None or getattr(ph,'integral',lambda:0)()<=0:
            rows.append([nm,lab,0.0,0.0]); block[nm][lab]={"photon_src_per_s_cm3":0.0,"contact_dose_uSv_h":0.0}; continue
        Es=np.asarray(ph.x); Sv=np.asarray(ph.p)          # E [eV], line strength [photons/s/cm3]
        src_tot=float(Sv.sum())
        mu=np.array([mu_of(mt,E) for E in Es])
        phi_s=Sv/(2.0*mu)                                  # surface flux photons/cm2/s
        h=np.interp(Es,de,dc)                              # pSv*cm2
        H_pSv_s=float(np.sum(phi_s*h)); H_uSv_h=H_pSv_s*3600e-6
        rows.append([nm,lab,src_tot,H_uSv_h])
        block[nm][lab]={"photon_src_per_s_cm3":src_tot,"contact_dose_uSv_h":H_uSv_h}
    print(f"{nm}: shutdown {block[nm]['shutdown']['contact_dose_uSv_h']:.3e} uSv/h  -> 1000yr {block[nm]['1000yr']['contact_dose_uSv_h']:.3e}",flush=True)

with open(os.path.join(RES,"shutdown_dose.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["region","cooling","photon_source_per_s_cm3","contact_dose_rate_uSv_per_h"])
    for r in rows: w.writerow([r[0],r[1],f"{r[2]:.4e}",f"{r[3]:.4e}"])

jp=os.path.join(RES,"designspec_results.json"); d=json.load(open(jp))
d["shutdown_dose"]={"model":"semi-infinite slab self-dose: phi_s(E)=S_v(E)/(2 mu(E)); H=sum phi_s*h_ICRP74_AP; units uSv/h; contact estimate not 3-D transport",
                    "per_region":block}
json.dump(d,open(jp,"w"),indent=2,default=float)
print("SHUTDOWN DOSE DONE",flush=True)
