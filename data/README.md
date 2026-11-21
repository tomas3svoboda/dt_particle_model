# Faner numerical records

This catalog contains **Faner cases only**. Superseded steam-rich PD/DT
recipes and their sensitivity/refinement results are withdrawn from main.
The current PD history is under `../revised_pd/`.

`final.json` preserves original Faner case values and three checksum-pinned
archive groups: `final_reference.tar.xz`, `final_adopted_faner.tar.xz` and
`final_faner_profile.tar.xz`. Downloads use the historical release and retain
its original archive/member hashes.

`ladder_report.json` is explicitly a Faner-only projection of the original
report. `variants/` retains Faner variants; `reproduction_verification.json`
preserves the completed Faner public-rerun check. No new qualification is
claimed by this cleanup.

```sh
python scripts/reproduce_final_reference.py --verify-manifest
python scripts/reproduce_final_reference.py --extract-only
```

The second command obtains only the three retained Faner archive groups.
Historical PD records remain versioned, but their case names are rejected by
the active public entry point.
