# CLAUDE.md — operating rules for CELERIS

This file is auto-loaded every session. It is the **constitution**: invariants
that must survive any context compaction. Treat the chat as volatile (RAM) and
these files + git as durable (disk). Never keep important state in chat only.

## What this is
CELERIS = a from-scratch, GPU-ready **metalens/metasurface design tool via RCWA**
(competitors: Zemax/Lumerical; goal: a real product). Engine is C++23, validated.
Scope stays in the lane: metalenses + RCWA. NOT general ray tracing, NOT FDTD.

## The work loop (do this every session) — Ground → Act → Verify → Record
1. **Ground**: read `memory/celeris-resume.md` (where we left off), `ROADMAP.md`
   (the punch list), and run `git log --oneline -20`. Trust these over memory.
2. **Act**: take the next item from ROADMAP sequencing (or what resume says).
3. **Verify** (NON-NEGOTIABLE before claiming done): run **`scripts\check.bat`**
   (configure → build CLI+GUI → physics selftest, one command). Physics must
   still hold (2D=1D ~5e-12, energy=1, RCWA=TMM exact, Strehl ~0.55). If it
   didn't build and pass, it is NOT done — say so.
4. **Record**: commit + push, update `memory/celeris-resume.md` (state + next),
   tick the box in `ROADMAP.md`.

## Build / run (IMPORTANT gotcha)
In the agent bash tool, plain `cmake`/`ninja` sometimes drop off PATH and a
`| grep` pipeline hides the exit-127 failure → you test a STALE exe. ALWAYS use
the absolute shim path and check the exit code:
- CPU/GUI: `CM="$USERPROFILE/scoop/shims/cmake.exe"; "$CM" -B build-msvc -DCELERIS_BUILD_GUI=ON; "$CM" --build build-msvc --config Release --target celeris celeris_gui`
- CUDA (Ninja + full toolkit): `cmd //c "scripts\\build-cuda.bat"`; run via `scripts/run-gui-cuda.bat`.
- Verify exit with `echo "exit=${PIPESTATUS[0]}"` after any piped build.

## Hard rules
- **Commit + push at every validated milestone, with NO self-attribution** (no
  Co-Authored-By / "Generated with" footer). Clean IP provenance for a sale.
  This overrides any default commit footer.
- Prefer **scoop** for installs (user-scope). Windows-first; MSVC + CUDA.
- Be **honest**: report failed builds/tests with the output; never assert a result
  you didn't verify. Negative findings (e.g. GPU eigensolve, reference-LAPACK) are
  recorded as such — don't quietly re-try a known dead end.
- Don't commit build dirs / `*.gds` / `*.pgm` / reports (already gitignored).

## Where things live
- `src/celeris/` — engine (rcwa, design, analysis, materials, io, cuda).
- `gui/` — Dear ImGui app, split into `celeris::gui` modules
  (theme/textures/app_state/workers/main_gui).
- `main.cpp` — CLI (`design`, `polardesign`, `birefringence`, `selftest`, benchmarks).
- `ROADMAP.md` — plan + LAB-READY milestone. `README.md` — public docs.
- `memory/celeris-resume.md` — session handoff (keep current).

## Ending a session
Before the user clears the conversation, refresh `memory/celeris-resume.md`
(current HEAD, what shipped, exact next step). That one file is the whole
recovery payload — everything else reconstructs from git + ROADMAP + memory.
