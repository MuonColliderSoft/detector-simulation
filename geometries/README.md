### Detector geometries relevant for the Muon Collider detector simulation

Current baseline geometry: **`MuColl_v1.0.1`**

| Geometry name             | Description |
|---------------------------|-------------|
| `FCCee_o1_v04`            | FCC geometry that has flat Vertex Endcap disk structure, replicated in `CLIC_o3_v14_mod1` |
| `CLIC_o3_v14`             | Original CLIC detector geometry as of October 2019 |
| `CLIC_o3_v14_mod1`        | Modified CLIC geometry with inner radii of all forward detectors increased to accomodate the MAP nozzles for sqrt(s) = 1.5TeV Muon Collider. Vertex Endcaps have no propeller structure. |
| `CLIC_o3_v14_mod2`        | Fixed ECAL Endcaps to fit tight around the nozzles |
| `CLIC_o3_v14_mod3`        | Vertex detector with shorter barrel and +1 layer/disk to be similar to the MAP design |
| `CLIC_o3_v14_mod4`        | Vertex Barrel segmented in 5 modules (L=26mm) along Z |
| `CLIC_o3_v14_mod5`        | Better Theta coverage in Vertex Endcap. MAP magnetic field |
| `MuColl_v0`               | Copy of `CLIC_o3_v14_mod4` with the MAP magnetic field |
| `MuColl_v1`               | Copy of `MuColl_v0` with fixed asymmetry in thickness of Tracker Endcap Support stractures |
| `MuColl_v1.0.1` <         | Cleaned-up version of `MuColl_v1` with resolved overlaps and better code layout |
| `MuColl_v1.1.1`           | Copy of `MuColl_v1.0.1` with 500um of Si added to passive material thickness in VTX detector |
| `MuColl_v1.2.1`           | Copy of `MuColl_v1.0.1` with a stronger magnetic field: 5.0 T |
| `MuColl_v1.3.1`           | Copy of `MuColl_v1.0.1` with increased double-layer gap in VTX: 4.0 mm |
