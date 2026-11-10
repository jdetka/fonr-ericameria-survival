SURVIVAL_RULE: last_alive

Interpretation:
  last_alive -> viability / establishment (alive through last confirmed alive; top-kill treated as temporary).

Primary inference choices (logic):
  - last_alive: Cox time-to-event (plant-level) + clustered SEs (patch),
               because it matches monotonic viability and avoids repeated-binary separation.
  - raw: GEE (plant-clustered) for marginal effects + cohort-interaction GEE,
         because raw includes top-kill + recovery and repeated 0/1 correlation is expected.

Standalone CAGE-only inference:
  - raw: CAGE-only GEE to isolate patch effects under herbivore exclusion.
  - last_alive: CAGE-only Cox to isolate patch effects on establishment without herbivory.

Contrasts (no emmeans):
  - Mechanism contrasts from full 3-way models:
      * Herbivory effect (NONE vs CAGE) within cohort × patch
      * Patch effect (OU vs IN) within cohort × treatment
  - Overall (treatment-averaged) ratios within each cohort:
      * Overall patch effect: geometric mean of patch effects under CAGE and NONE
      * Overall herbivory effect: geometric mean of herbivory effects within IN and OU

Supporting models:
  - Scaled GLMM (may warn under separation): used for comparison/sensitivity.
  - GLMM without plant RE: diagnostic stability check (often fits faster/cleaner).
  - Site-fixed GLM: checks that treatment/context patterns persist controlling for site.
  - Optional cohort-interaction GLMM (NO plant RE): only if explicitly enabled.

Run generated: 2026-02-23 17:14:04.491875
