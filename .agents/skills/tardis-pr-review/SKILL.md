---
name: tardis-pr-review
description: Analyze a numbered pull request in tardis-sn/tardis and produce a concise, plain-language guide to changed equations, physical assumptions, scientific data, numerical methods, and supporting tests that merit close review by scientists. Use when the user supplies a TARDIS PR number and asks for scientific review priorities or areas of scientific interest. Ask only questions directly motivated by changed code or its immediate scientific effects, and explain their importance without assuming Git or software-review expertise. Do not use for a comprehensive code review, implementation work, local-only changes, non-TARDIS repositories, or requests without a PR number.
---

# Highlight Science Questions in a TARDIS Pull Request

## Establish a Read-Only Snapshot

1. Require a positive integer pull-request number. Treat `tardis-sn/tardis` as
   the canonical repository and work from a TARDIS checkout.

2. Record the current branch, HEAD, remotes, and `git status --short`. Preserve
   all user work. Do not switch branches, stash, reset, clean, pull, or use
   `gh pr checkout`.

3. Obtain the PR title, body, URL, state, author, base ref and SHA, head SHA,
   labels, commits, and changed files. Prefer:

   ```shell
   gh pr view <number> --repo tardis-sn/tardis \
     --json number,title,body,url,state,isDraft,author,baseRefName,baseRefOid,headRefName,headRefOid,labels,commits,files
   ```

   Read linked requirements and relevant existing discussion. When using the
   REST API, paginate commits, files, reviews, review comments, and issue
   comments instead of trusting first-page results. Do not install missing
   tools.

4. Fetch the base and head into review-only refs:

   ```shell
   git fetch --no-tags https://github.com/tardis-sn/tardis.git \
     +refs/heads/<base-ref>:refs/tardis-pr-review/<number>/base \
     +refs/pull/<number>/head:refs/tardis-pr-review/<number>/head
   ```

   Verify the fetched object IDs against the recorded SHAs. Re-query once if
   the base moved; report an inconsistent snapshot if either SHA still differs.
   Never execute code from the fetched PR.

## Recover the Scientific Intent

1. Follow the active root `AGENTS.md`. Read applicable current developer
   documents from the fetched base tree, especially feature workflow, design
   philosophy, and testing and regression guidance. Treat them as review policy,
   not evidence that the proposed science is correct.

2. Read the PR requirement, complete triple-dot diff, and relevant code, tests,
   fixtures, configuration, physics documentation, cited equations, and data
   provenance from the base and head objects. Prefer `git diff`, `git show`, and
   `git grep` with external diff and text-conversion drivers disabled so the
   working tree remains unchanged.

3. State the claimed scientific outcome, the relevant physical regime, and the
   explicit and implicit assumptions. Trace the changed quantities through
   their earlier inputs and later uses—for example, from configuration, model,
   or atomic data through plasma, opacity, transport, estimators, convergence,
   and spectrum generation. Do not limit the trace to changed files.

## Select Areas That Need Scientific Attention

Include an area only when a change, or an assumption connecting calculations,
can alter a scientific result, its interpretation, or the evidence used to
validate it. Scan the entire diff, but do not turn the brief into a file-by-file
software review.

Require this evidence chain for every reported area:

1. An exact changed line.
2. The scientific quantity, assumption, or test affected by that line.
3. A plausible effect on a calculation or result.
4. A specific unanswered question that would help a scientific reviewer assess
   that effect, or a directly established scientific consequence.

Omit the area if any link is missing. Use the following lenses only to interpret
actual changes; never turn them into a generic checklist:

- Physical formulation: equations, conservation laws, approximations, regimes
  of validity, constants, dimensions, units, and normalization.
- Frames and geometry: comoving versus observer quantities, frequency and
  wavelength transformations, time since explosion, homologous-expansion
  assumptions, shell boundaries, and radial or directional conventions.
- Numerical representation: discretization, interpolation, integration,
  tolerances, convergence, stability, floating-point behavior, array shapes,
  masks, ordering, indexing, and boundary handling.
- Monte Carlo transport: packet initialization and propagation, interaction
  probabilities, random sampling and seeds, estimators, packet weights, energy
  accounting, statistical uncertainty, and reproducibility.
- Matter-radiation coupling: ionization and excitation populations, LTE/NLTE
  assumptions, radiation-field updates, opacities, rates, Sobolev and macro-atom
  treatment, continuum processes, iteration ordering, and cached state.
- Scientific data semantics: atomic or decay-data provenance, identifiers,
  level and transition mappings, units, missing values, filters, joins, and
  conversions between tabular and array representations.
- Observable construction: spectral binning, luminosity integration, virtual or
  formal-integral spectra, escaping-packet selection, and output quantity
  definitions.
- Scientific validation: analytic or limiting cases, independent comparisons,
  regression baselines, tolerances, packet counts, stochastic variance, test
  oracle provenance, and assertions on physically meaningful quantities.

Assign each area exactly one of these labels, choosing the most specific label
supported by the changed lines:

- `Physical model change`: equations, physical assumptions, conservation laws,
  approximations, or regimes of validity change.
- `Convention change`: a frame, geometry, sign, normalization, unit, indexing,
  or other scientific convention changes.
- `Coupled calculation with downstream effects`: a change in one calculation
  changes inputs, state, or iteration behavior in later coupled calculations.
- `Change in interpretation of scientific data`: the meaning, provenance,
  mapping, filtering, units, or conversion of scientific data changes.
- `Change in stored reference result`: an expected or regression result changes,
  whether because the calculation or the stored comparison changes.
- `Calculation details changed`: discretization, interpolation, tolerances,
  ordering, floating-point handling, or other numerical details change without
  a more specific close-review classification above.
- `Defaults changed`: a default parameter, option, packet count, seed, or other
  implicit user input changes.
- `Boundary values changed`: shell, grid, mask, cutoff, range, or other
  boundary handling changes.
- `Tests changed`: tests, fixtures, assertions, tolerances, or test data change
  in a way that can alter or conceal scientific validation.

Use the first five labels for changes needing close review and the last four
for changes to check carefully. If a change genuinely covers multiple distinct
labels, split it into separate areas. These labels show where scientific
attention is most useful; they do not mean that the change is wrong.

For each area:

1. Link to the relevant changed lines in the version reviewed and name important
   later calculations when needed.
2. Explain the changed scientific quantity or assumption and the results it can
   affect later in the simulation.
3. Identify the physical regime, relation that should remain true, convention,
   or meaning of the data that a reviewer should verify.
4. Ask only questions caused by the changed lines or their immediate effects.
   Make each question identify the changed choice and the scientific uncertainty
   it creates. Do not ask broad questions such as whether units, conservation,
   tests, or documentation are generally adequate.
5. Distinguish demonstrated evidence from inference. If an error is directly
   established, describe its scientific consequence without inventing a
   question; otherwise present the uncertainty as a review question, not a
   defect.
6. Prefer one question and never include more than three questions for one area.
   Combine overlapping areas and remove duplicate questions.
7. Do not ask a question already answered by the changed code, tests, PR
   description, or review discussion. State the answer briefly only when it
   materially affects the remaining scientific question.

Exclude style, naming, typing, general maintainability, documentation polish,
and unrelated CI failures unless they obscure scientific meaning or validation.
Do not assign software-bug severity labels, give an overall approval verdict, or
imply that passing tests establish scientific correctness.

## Inspect Relevant Validation Evidence

1. Do not run PR code, tests, Ruff, Python imports, notebooks, documentation
   builds, or benchmarks locally. Inspect exact-head CI and existing artifacts
   only for evidence tied to a reported area.

2. Determine whether relevant unit, integration, regression, GPU, or benchmark
   jobs ran against the recorded head or its current synthetic merge commit.
   Treat skipped, label-gated, stale, cancelled, or absent checks as missing
   evidence. Inspect failed logs when they clarify a scientific validation gap.
   Do not summarize unrelated checks.

3. Never post a review or comment, approve, label, rerun, or otherwise mutate
   GitHub state without separate user authorization.

4. Re-query the PR head SHA before reporting. If it changed, fetch the new
   snapshot and inspect the incremental diff before relying on prior analysis or
   CI.

## Write for a Scientific Reader

Assume the reader understands the relevant astrophysics but may not know Git,
GitHub Actions, or software-review terminology.

1. Use precise scientific language and ordinary prose. Explain abbreviations on
   first use unless they are standard in the immediate TARDIS context.

2. Translate software terms in the report:

   - Say **version reviewed**, not `head`, `SHA`, or `commit`, except when giving
     the identifier in the final version details.
   - Say **automated tests**, not `CI`, `check`, or `job`.
   - Say **stored reference result** or **expected result**, not `regression
     baseline` or `test oracle`.
   - Name the actual later calculation or scientific output instead of saying
     `caller`, `consumer`, or `downstream`.
   - Describe the expected inputs, outputs, units, and meanings instead of using
     `interface` or `data contract`.
   - Say **area needing scientific attention**, not `hotspot`.
   - Introduce a source link as **Where to look**, not a `line anchor`.

3. When a technical software term is unavoidable, define it in one short
   sentence at first use. Do not expect the reader to interpret branch names,
   workflow names, test markers, or abbreviated check names.

4. Focus each area on the scientific choice and its possible effect. Mention
   code mechanics only when they help the reader understand that effect.

5. Keep the brief short. Do not repeat the PR description, explain familiar
   TARDIS physics, narrate the inspection process, list unaffected files, or
   restate the same consequence under several headings.

## Report a Scientific-Review Brief

Lead with: “This brief highlights science questions for specialist review; it
does not determine whether the proposed change is correct.” Follow it with at
most one sentence describing what the PR aims to change.

Apart from that opening and the final version line, report only **Areas needing
scientific attention**. Give each area a scientific title and exactly one of the
specific labels above, followed by:

- **Where to look**: link to the relevant lines.
- **Why it matters**: explain the possible scientific effect in one or two
  sentences.
- **Question for the reviewer**: include this only for an unresolved question.
  Ask one concise question directly tied to the changed lines. Add a second or
  third only when each addresses a distinct scientific uncertainty created by
  those changes. Do not invent a question for an established consequence.
- **Evidence**: include at most one sentence when an automated test, comparison,
  or missing result changes how the reviewer should approach the question.

Add one short paragraph about how effects move through the simulation only when
it is necessary to connect several reported areas. Do not include a separate
general test summary.

End with one **Version reviewed** line containing the PR URL and identifiers for
the starting and proposed versions, followed by a note that the proposed code
was inspected but not run locally.

If no changed area merits scientific review, say so explicitly and briefly
explain why the change appears science-neutral, then give the version reviewed.
Do not fill the report with general software-review observations.
