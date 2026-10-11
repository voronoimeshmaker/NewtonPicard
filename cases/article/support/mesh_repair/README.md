# Local VMM exact-clipping repair — 2026-10-10

The installed VoronoiMeshMaker rejects the random meshes with seed 12345 on
[-15,15] x [-5,5]: L3 (87552 sites, site 902 edge 3) and L4 (349184 sites,
site 877 edge 0). For L3 the nearest reversed face exists in cell 20454;
endpoint mismatch is 7.040162925441755e-13 versus matching tolerance
2.1316282072803006e-13. Original face length is 0.055845280869207432.

## Repair scope

MeshFactory.cpp is a local replacement for the historical MeshFactory translation
unit. It first uses the unchanged installed builder. Only the identified
reciprocal-segment exception triggers a complete exact reconstruction from the
same sites and boundary. All other exceptions propagate. The installed library
and historical copied sources are unchanged.

ArticleExactVoronoiCellBuilder2D computes bisectors and intermediate polygon
intersections with CGAL exact constructions from the original binary floating
point site coordinates. It forces exact evaluation before converting final
vertices to double. It retains the original duplicate-vertex policy. The renamed
clipped-diagram builder uses this cell builder and retains the installed
connectivity checks, including tolerance, unique reciprocal partners, face length
and geometry validation. No face is silently deleted and no matching threshold
is relaxed. Solver settings, random seeds and sample sets are unchanged.

## Verification

Both full random meshes were accepted by the original connectivity checks:
- L3: 87552 cells, 522938 directed internal faces, total-area error 2.50e-15.
- L4: 349184 cells, 2090150 directed internal faces, total-area error 6.11e-16.

The small-mesh regression test compares original/exact neighbour sets and cell
areas and deliberately moves a shared vertex to confirm broken connectivity is
still rejected. Source: test_exact_clipping.cpp.

The 200 Cartesian/hexagonal results at each refinement remain in their original
run directories. Resumed configurations contain only random. pending_cases.py
skips already completed configured campaigns; scientific nonconvergence remains
a recorded outcome, distinct from a setup failure.

## Follow-up in VoronoiMeshMaker

Port robust construction into the owning source repository, with regression
coverage for both seeds/resolutions, boundary and near-degenerate cases, and
performance measurements. Do not merely enlarge the face-pairing tolerance.
This local repair is not an installed VMM release or a full VMM-wide audit.
