# CLAUDE.md — operating rules for CELERIS

This file is auto-loaded every session. It is the **constitution**: invariants
that must survive any context compaction. Treat the chat as volatile (RAM) and
these files + git as durable (disk). Never keep important state in chat only.

## What this is
CELERIS = a from-scratch, GPU-ready **metalens/metasurface design tool via RCWA**.
Engine is C++23, validated. Scope stays in the lane: metalenses + RCWA. NOT general
ray tracing, NOT FDTD.

**Primary goal (current, changed 2026-10-19): publish the CELERIS code paper in the
_Journal of Open Source Software_ (JOSS)**, solo (undergrad author), + an arXiv
long-form report. CPC was dropped for the code paper (novelty bar too high for a solo
first paper); SciPost bars solo undergrads. `paper/paper.md` is the JOSS paper;
`paper/longform/` is the arXiv companion (ex-CPC manuscript), NOT a second journal paper.
**Earliest JOSS submission: 2026-12-17** (repo public 2026-06-16 + >6 months). The real
blockers are JOSS's desk-reject gates — sustained (not bursty) development, public
issues/PRs, concrete research impact, API docs — tracked in `ROADMAP.md` §J. Commit real
work steadily, through issues + PRs; **never fake or backdate history.** Later papers
(method/adjoint, results) are follow-ons — **one paper per distinct content, never the
same paper to multiple venues.**
**Read `memory/celeris-resume.md` for the full handoff** (venue rationale, session-scratchpad
warning, exact next steps). The Zemax/Lumerical *product* is the "PARKED" backlog in
`ROADMAP.md`, not the driver. CELERIS is **Apache-2.0** (permissive; keeps commercial
optionality).

## The work loop (do this every session) — Ground → Act → Verify → Record
1. **Ground**: read `memory/celeris-resume.md` (where we left off), `ROADMAP.md`
   (the punch list), and run `git log --oneline -20`. Trust these over memory.
2. **Act**: take the next item from ROADMAP sequencing (or what resume says).
3. **Verify** (NON-NEGOTIABLE before claiming done): on Linux run
   `cmake -B build-linux -G "Unix Makefiles" -DCELERIS_BUILD_GUI=OFF && cmake --build
   build-linux --target celeris -j16 && ./build-linux/celeris selftest` (exit 0; on
   Windows: **`scripts\check.bat`**, configure → build CLI+GUI → physics selftest). Physics must
   still hold (2D=1D ~5e-12, energy=1, RCWA=TMM exact, Strehl ~0.55). If it
   didn't build and pass, it is NOT done — say so.
4. **Record**: commit on the branch → PR → CI green → merge (see JOSS process
   rules); update `memory/celeris-resume.md` (state + next), tick the box in `ROADMAP.md`.

## JOSS process rules (EVERY session — these are the submission gates)
JOSS desk-rejects on *process*, not just code (ROADMAP §J). Each session must leave
the repo looking like a steadily maintained open project, because it is one.
- **Steady development, never bursts.** One issue (or one coherent slice of one) per
  session, a few focused commits. Do not dump a week of work into one day; if a task
  is big, split it into PRs on different days. **Never backdate, fake, or pad commits**
  (no empty/"activity" commits) — real work only.
- **Everything goes through issues + PRs.** Start from a GitHub issue (open one if
  none fits), work on a branch, open a PR whose body says `Closes #N`, wait for CI
  green, then merge (`gh pr merge --merge --delete-branch`). No direct pushes to `main`
  except handoff/memory updates. Label issues (`enhancement`, `bug`, `documentation`).
- **API docs are part of "done".** Every public declaration in
  `src/celeris/**/*.hpp` gets a Doxygen `///` comment (`///<` for members); every
  pybind11 `m.def`/`py::class_` gets a docstring. `doxygen docs/Doxyfile` must stay
  warning-free (CI enforces it). A change that adds an undocumented public API is
  not done.
- **Research impact must be real and public.** As of 2026-10-19 there is no lab/PI
  use — the author uses CELERIS on their own. The plan is to *create* honest evidence:
  (1) a real CELERIS result posted as a preprint (the single-etch group-delay-limit
  study, see ROADMAP R1), (2) a poster/talk (e.g. UCF's undergrad research showcase),
  The project is SOLO by design — do not suggest collaborators, testers, or outreach.
  Log each concrete use (date, project, what CELERIS computed, any talk/poster/
  preprint) in `memory/research-use.md` **only from what the author tells you — never
  invent uses**. That log feeds the paper's Research impact statement. The repo is
  public: record only what the author says may be public.

## Build / run (IMPORTANT gotcha)
In the agent bash tool, plain `cmake`/`ninja` sometimes drop off PATH and a
`| grep` pipeline hides the exit-127 failure → you test a STALE exe. ALWAYS use
the absolute shim path and check the exit code:
- CPU/GUI: `CM="$USERPROFILE/scoop/shims/cmake.exe"; "$CM" -B build-msvc -DCELERIS_BUILD_GUI=ON; "$CM" --build build-msvc --config Release --target celeris celeris_gui`
- CUDA (Ninja + full toolkit): `cmd //c "scripts\\build-cuda.bat"`; run via `scripts/run-gui-cuda.bat`.
- Verify exit with `echo "exit=${PIPESTATUS[0]}"` after any piped build.

## Hard rules
- **Commit + push at every validated milestone, with NO self-attribution** (no
  Co-Authored-By / "Generated with" footer). Clean, single-author provenance — for
  the papers (JOSS/CPC authorship) and any future sale. Overrides any default footer.
- Prefer **scoop** for installs (user-scope). Dev is Windows-first (MSVC + CUDA), but
  **a clean Linux/CPU-only build + CI is now top priority** — JOSS/CPC reviewers run
  it on Linux without a GPU (ROADMAP §1).
- Be **honest**: report failed builds/tests with the output; never assert a result
  you didn't verify. Negative findings (e.g. GPU eigensolve, reference-LAPACK) are
  recorded as such — don't quietly re-try a known dead end.
- Don't commit build dirs / `*.gds` / `*.pgm` / reports (already gitignored).

## Where things live
- `src/celeris/` — engine (rcwa, design, analysis, materials, io, cuda).
- `gui/` — Dear ImGui app, split into `celeris::gui` modules
  (theme/textures/app_state/workers/main_gui).
- `main.cpp` — CLI (`design`, `polardesign`, `birefringence`, `selftest`, benchmarks).
- `ROADMAP.md` — publication plan (JOSS/CPC milestones §0–§3) + PARKED product
  backlog. `README.md` — public docs.
- `memory/celeris-resume.md` — session handoff (keep current).

## Ending a session
Before the user clears the conversation, refresh `memory/celeris-resume.md`
(current HEAD, what shipped, exact next step). That one file is the whole
recovery payload — everything else reconstructs from git + ROADMAP + memory.
