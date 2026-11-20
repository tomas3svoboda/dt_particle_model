"""Resolved single-particle model (Phase 1) — the frozen GT-PS-2 sphere.

Built in sub-gate order (see ``docs/CORE2_PHASE1_PARTICLE.md``), each mechanism
verified in isolation with a manufactured-solution order check before coupling,
to avoid over/under-determining the equation set:

    grid        spherical finite-volume geometry
    dry_shell   PHY-019 Fick pore-vapor + PHY-020 local-equilibrium storage   [gate 1c]
    wet_core    PHY-013 conduction + PHY-038 shell-conserved n-hexane          [gate 1d]
    front       PHY-039 Rankine-Hugoniot mass+energy jumps + fugacity          [gate 1e]
    dry_thermo  dry-region binary storage + common-datum energy                [gate 1f]
    dry_region  conservative annular ALE component/energy ledgers              [gate 1f]
    sphere      high-loading state + exact lossless activation event            [gate 1f]
    sphere_coupling  coupled conservation/event-consumption qualifiers          [gate 1f]
    thermal_oracle  fixed-sphere finite-film radial/lumped qualification         [gate 1g]
    coupled_pore    G1G-01 stationary-water / binary-MS DAE foundation           [gate 1g]
    conditioning    unit-invariant physical-Jacobian rank/condition certificate   [gate 1g]
    bed_film        correlated nonqualifying external heat/mass film oracle       [gate 1g]
    transport_coefficients  bounded nonqualifying binary face mobility           [gate 1g]
    water_topology_oracle  stationary-versus-front-fitted manufactured audit      [gate 1g]
    coupled_transport  fully dry all-species/energy predictive solve              [gate 1g]
    cut_geometry    fixed material grid with exact spherical cut cells            [gate 1g]
    wet_water       wet retained-water/energy predictive solve                    [gate 1g]
    cut_transport   square fixed-cut ALE/RH residual foundation                   [gate 1g]
    cut_integrator  safeguarded conservative same-cell nonlinear solve             [gate 1g]
    cut_birth_event / cut_birth_integrator  exact A1 dry-shell birth               [gate 1g]
    cut_face_event / cut_face_event_integrator  exact master-face arrival          [gate 1g]
    cut_face_to_face  exact adjacent-master-face unknown-duration residual          [gate 1g]
    cut_face_departure / cut_face_departure_integrator  fail-closed seam audit      [gate 1g]
    cut_face_tangent  exact zero-volume departure tangent/overlap audit             [gate 1g]
    cut_event_orchestrator  atomic one-interior-face requested macrostep             [gate 1g]
    cut_extinction_event / cut_extinction_projection  exact supplied-time endpoint [gate 1g]
    cut_extinction_time  positive-core finite-time localization certificate        [gate 1g]
    cut_extinction_orchestrator  atomic exact-dry remainder transfer                [gate 1g]

Structural no-contradiction invariants (enforced by module boundaries):
  - the wet core carries NO diffusion operator (shell-conserved, PHY-038);
  - the dry shell carries the Fick operator (PHY-019); the two never overlap;
  - the front position s comes from the RH jumps ONLY, never Faner s(X);
  - latent/binding energy enters ONLY through common-datum state differences,
    never as a separate j*dHvap source (FLUX-ENERGY-INTERFACE).
"""

from __future__ import annotations

__all__ = [
    "bed_film",
    "conditioning",
    "coupled_pore",
    "coupled_transport",
    "cut_birth_event",
    "cut_birth_integrator",
    "cut_continuation",
    "cut_event_orchestrator",
    "cut_extinction_event",
    "cut_extinction_orchestrator",
    "cut_extinction_projection",
    "cut_extinction_time",
    "cut_face_departure",
    "cut_face_departure_integrator",
    "cut_face_event",
    "cut_face_event_integrator",
    "cut_face_tangent",
    "cut_face_to_face",
    "cut_face_to_face_integrator",
    "cut_geometry",
    "cut_integrator",
    "cut_sparsity",
    "cut_transport",
    "dry_region",
    "dry_shell",
    "dry_thermo",
    "front",
    "grid",
    "local_equilibrium_diagnostic",
    "sphere",
    "sphere_coupling",
    "surface_active_set",
    "thermal_oracle",
    "transport_coefficients",
    "water_topology_oracle",
    "wet_core",
    "wet_water",
]
