# FvFramework4Students (NTFD Edition)

FvFramework4Students v1.1c is a MATLAB teaching framework for the course Numerical Techniques in Fluid Dynamics (NTFD). It provides the scaffolding to implement cell‑centred finite volume (FV) discretisations in 2‑D: mesh/domain geometry, field containers, a sparse linear‑system bridge, and plotting utilities. You write the physics (coefficients, boundary conditions, iteration); the framework handles the bookkeeping.

This README expands on the framework in the context of the NTFD sessions and the accompanying slides (Session 1: Practicalities & MATLAB framework; Session 2: Steady‑state Diffusion).
GOOFY MFER

## Quick Start

- Launch MATLAB in this folder or double‑click `launch.bat` (runs `launchfvframework4students`).
- MATLAB prints the license notice, adds `src` and subfolders to the path, then moves into the `work` directory.
- Run `examplecase/runexamplecase` to generate a simple mesh, set fields and execute the example solver. Use it as a template for your own cases.

Target environment: MATLAB with `classdef` support (R2014+). No extra toolboxes required.

## Course Context (from the Session Slides)

- Prerequisites: fluid mechanics, numerical modeling, MATLAB.
- Format: primarily self‑study with weekly contact sessions (blackboard theory/explanations).
- Assessment: 70% project work, 30% presentation to peers.
- Final report (deliverables):
  - Describe discretisation choices (spatial and temporal, if applicable).
  - Perform grid convergence and order‑of‑accuracy estimation.
  - Verify results (analytically or against reference solutions) and discuss.
  - Keep it concise and to the point. Submit the PDF to the instructors (timing around start of January exam period; confirm exact date in class).

### Suggested Project Milestones (checkpoints)

1) Geometry implementation of normals, tangents, etc.
2) 2‑D steady‑state diffusion.
3) 2‑D steady‑state convection–diffusion.
4) 2‑D channel flow with imposed pressure field.
5) SIMPLE pressure correction for 2‑D channel flow.
6) Rhie–Chow interpolation to avoid checkerboarding.
7) 2‑D lid‑driven cavity.
8) … extend as needed for your project.

These align naturally with the framework components described below.

## Repository Layout

- `launchfvframework4students.m`, `launch.bat` – entry points: print the licence, add `src` to the MATLAB path, and enter `work`.
- `src/auxiliary` – infrastructure: `FvmLabEntity`, `FvmLabEntityList`, `TypeChecker` validation helpers, singly linked list utilities, `StopfileMonitor`.
- `src/domain` – geometry and topology: `FvMesh` (raw mesh data and validation), `FvDomain` (derived metrics and default zones), zone classes (`@CellZone`, `@FaceZone`, `@BfaceZone`, `@RangeZone`, `@IZone`, `@VertZone`), and numerical helpers in `domain/num` (connectivity maps, volumes, normals, ghost‑cell mirroring, etc.).
- `src/fields` – `Value` (algebra/storage for scalar/vector/tensor data), `Field` (binds a `Value` to a zone), factories `newfield.m` / `newvalue.m`, plus low‑level kernels in `L1`, `L2` for arithmetic and dot products.
- `src/linear` – `ScalarFvEqn2`: an FV‑oriented sparse system container that maps stencil coefficients to MATLAB `sparse` matrices.
- `src/meshing` – `LineSeed` (biased 1‑D distributions) and `TwoSeedMesher.genmesh` (orthogonal block meshes with orderly numbering and named boundaries).
- `src/solvers` – example skeleton(s), e.g. `examplesolver.m` showing a residual‑controlled loop and conversion to `sparse`.
- `src/visualize` – `fvmplotmesh`, `fvmplotfield`, vector overlays, and index labelling for inspection.
- `src/tools/install` – path helpers used by the launcher (`addtopath`, `fvgenpath`).
- `work/` – your sandbox: ships with `examplecase/runexamplecase.m` wiring mesh → domain → fields → solver → plots.

## Coefficient Storage Layout (ScalarFvEqn2)

`ScalarFvEqn2` stores the linear system A·x = b in a compact 1‑vector `adata` that follows the mesh topology. The ordering is:

- First the `n` diagonal entries A(i,i) for i = 1..n.
- Then, for every internal face j = 1..nIf:
  - A(fn1(j), fn2(j)) and A(fn2(j), fn1(j)) occupy the next two entries in `adata`.
- Then, for every boundary face j = nIf+1..nF (coupling owner cell to its ghost cell):
  - The pair A(fn1(j), fn2(j)) and A(fn2(j), fn1(j)) occupy the final 2·nBf entries in `adata`.

The right‑hand side is the dense vector `bdata` of length `n`. Use `[A, b] = to_msparse(eqn)` to materialise MATLAB’s `sparse` matrix.

## Data Flow and Core Workflow

1) Mesh generation – Create 1‑D seeds with `LineSeed.lineSeedOneWayBias` and call `TwoSeedMesher.genmesh(seedI, seedJ, boundaryNames)` to form a structured block. Boundaries receive names (e.g. `WESTRAND`, `OOSTRAND`, `ZUIDRAND`, `NOORDRAND`).
2) Domain creation – Build a `FvDomain` via `newdomain(mesh, 'MyDomain')`. The domain computes:
   - Connectivity (vertex/face/cell neighbours and ordered cell faces).
   - Geometry (centroids, normals, areas, volumes, face–cell distances).
   - Ghost cells for boundary faces by mirroring owner‑cell centroids.
   - Default zones: `allCells`, `allFaces`, `allBfaces`, `allVerts`, plus any named boundary zones from the mesh.
3) Fields – Construct `Field(zone, rank)` on zones (rank 0 = scalar, 1 = vector, 2 = tensor) and set/reset their values. Algebra on fields is supported through low‑level kernels to keep operations consistent.
4) Equation assembly – Use `ScalarFvEqn2(dom)` to hold a system in compressed, topology‑aware form. Fill:
   - `adata`: diagonal terms for each cell, then off‑diagonals per face (internal first, then boundary; both neighbour directions).
   - `bdata`: right‑hand side.
   Convert to MATLAB with `[A, b] = to_msparse(eqn)`.
5) Solve & iterate – Solve with `A \ b` (direct) or Krylov methods; compute residuals; update fields; check tolerances and iteration caps.
6) Post‑process – Use `fvmplotmesh`, `fvmplotfield`, and restriction to zones for diagnostics. Label indices when validating stencils.

## MATLAB Workflow Tips

- Use `clear variables` to reset workspace variables without flushing compiled classes.
- Use `clear classes` after editing class definitions to force MATLAB to reload.
- `close all`, `clc` help keep sessions tidy.
- In the MATLAB editor, press Ctrl+D while the cursor is on a function call to open that function quickly.

## Session 2 Focus: Steady‑State Diffusion

The transport equation for a passive scalar ϕ is (generic form):

- d(ϕ)/dt + div(u ϕ) = div(Γ grad ϕ) + S

Session 2 simplifies to the steady diffusion problem:

- Steady state (d(ϕ)/dt = 0)
- No advection (div(u ϕ) = 0)
- Uniform diffusivity Γ [m^2 s^-1]
- Governing equation reduces to: div(Γ grad ϕ) = S

Implementation guide in this framework:

- Create `eqn = ScalarFvEqn2(dom)`.
- For each internal face j between cells P = fn1(j) and N = fn2(j):
  - Compute diffusive coupling a_diff,j from Γ, face area and cell‑to‑cell distance (use domain geometry such as `fArea`, `fXi`, `fXiMag`, `fDs`).
  - Accumulate to `adata` so that A(P,N) and A(N,P) get −a_diff,j and the diagonals A(P,P), A(N,N) get +a_diff,j.
- For boundary faces, use the domain’s ghost‑cell convention (each boundary face couples the owner cell to its ghost cell) to impose boundary conditions by appropriate diagonal/off‑diagonal and RHS contributions.
- After filling `adata` and `bdata`, call `[A, b] = to_msparse(eqn)`, solve, and write back into a `Field` with `set`.

Tip: The provided `work/examplecase/runexamplecase.m` already sets up a mesh and Dirichlet‑like data. Use it to validate your diffusion implementation and plotting.

## Boundary Conditions and Ghost Cells

- The domain adds a mirrored ghost cell for each boundary face. Boundary faces connect (c1) to the physical owner cell and (c2) to a ghost cell indexed after the physical cells. This supports popular FV strategies for Dirichlet/Neumann enforcement.
- Named boundary zones from meshing (e.g. `WESTRAND`) are accessible via `getzone(dom, 'WESTRAND')` and can be used to assemble boundary contributions selectively.

## Example: Minimal Diffusion Solver Skeleton

Use `src/solvers/examplesolver.m` as a starting point. The essential steps are:

```matlab
function result = mydiffusionsolver(casedef)
dom = casedef.dom;
phi = Field(dom.allCells,0);  % scalar field
reset(phi, 0);

eqn = ScalarFvEqn2(dom);
iterate = true; niter = 0; 
while iterate
  niter = niter + 1;
  reset(eqn);
  % 1) Loop faces/cells: assemble diffusion contributions into eqn.adata
  % 2) Apply sources into eqn.bdata
  % 3) Apply BCs via ghost-cell rows/cols and eqn.bdata
  [A,b] = to_msparse(eqn);
  x = A\b; set(phi, x');
  % 4) Check residual and stop criteria
  res = norm(b - A*x);
  if res < casedef.iteration.TTol || niter >= casedef.iteration.maxniter
    iterate = false;
  end
end
result.phi = phi;
end
```

Focus on correctness of face‑based couplings and consistent treatment of boundary faces.

## Visualisation Tips

Common patterns:

```matlab
figure; hold on; axis equal; colormap(jet(50));
fvmplotmesh(casedef.dom, 1);
fvmplotfield(result.T, 'lin', 1);            % scalar contours
% Restrict a vector field to a boundary zone and plot arrows:
Uoost = restrictto(U, getzone(casedef.dom, 'OOSTRAND'));
fvmplotvectorfield(Uoost, 1);
```

When debugging, label indices with `fvmplotcellnumbers`, `fvmplotfacenumbers`, `fvmplotvertexnumbers`.

## Key Files to Read First

- `src/domain/@FvMesh/FvMesh.m` – minimal mesh representation and validation (`TypeChecker` enforced).
- `src/domain/@FvDomain/FvDomain.m` – topology, geometry and zone construction; ghost‑cell handling.
- `src/fields/@Field/Field.m` and `src/fields/@Value/Value.m` – field storage, algebra and restriction.
- `src/linear/@ScalarFvEqn2/ScalarFvEqn2.m` – how `adata`/`bdata` map to `sparse` A, b.
- `src/meshing/@TwoSeedMesher/genmesh.m` – numbering conventions and named boundary ranges.
- `work/examplecase/runexamplecase.m` – end‑to‑end wiring and plotting.

## Extending the Framework

- New cases – Create a subfolder under `work/`, copy `runexamplecase.m` and iterate.
- Custom solvers – Start from `examplesolver.m`; adopt a clear assembly loop and residual checks.
- Boundary conditions – Use named zones to target boundaries; enforce via ghost‑cell couplings and RHS terms.
- Meshes – Extend `TwoSeedMesher` or produce structs compatible with `FvMesh.fromstruct`, respecting ordering (internal first).
- Utilities – Reuse `TypeChecker` for input validation and `FvmLabEntityList` for collections.

## Course Deliverables Checklist

- Discretisation rationale stated (and linked to equations used).
- Grid‑refinement study with order‑of‑accuracy estimate.
- Verification against analytical or trusted numerical references.
- Clear figures (mesh, fields, convergence). Concise narrative.

## Troubleshooting

- After editing class interfaces, run `clear classes` in MATLAB.
- If paths aren’t found, re‑run `launchfvframework4students` to refresh path setup.
- Use index plotting utilities to diagnose stencil assembly issues.

## Licence

`license.txt` restricts usage to K.U.Leuven students and personnel for didactical purposes and personal study. Ensure compliance before redistribution or integration elsewhere.
