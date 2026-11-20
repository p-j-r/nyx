# nyx market note (2026-11-20)

Desk research from web searches done on 2026-11-20. It answers "what is the biggest use case for nyx, and who else is there". Nothing here is customer research: no one who would pay has been asked. Each claim carries a confidence label: **sourced** (a page says it), **inferred** (my reasoning from sourced facts), **unchecked** (I did not verify it).

## 1. Short answer

1. **The AI-agent design loop for photonics is not an empty field.** Flexcompute (Tidy3D, PhotonForge) already ships an MCP assistant (Nov 2025), published an autonomous design loop with Claude or Codex as the agent (June 2026), released PhotonForge 2.0 "built for engineers and AI agents" (Oct 2026) and took strategic investment from Samsung Ventures and Applied Ventures (Oct 2026). COMSOL has published agentic-workflow guides. "An AI can drive a solver" is no longer a differentiator.
2. **The "coupled optical and mechanical" gap is thinner than I first wrote.** The gap check (section 4) found that an agent already composes it from separate tools: on 17 Nov 2026 Flexcompute published an acousto-optic modulator design where an AI agent (Opus 5.5) used Tidy3D for the optics, the open-source NGSolve library for the piezoelectric acoustics and PhotonForge for layout. COMSOL does it in one tool today. What nyx offers is one local tool with one mesh, one set of units and a validation suite, not a capability nobody has.
3. **The market that pays is industry PIC design, the market that cites is academia.** Academic licences are cheap, so the free price is not the pitch there. The pitch is agent-readiness and trust.
4. **The thing worth building and publishing is trust**: the validation suite. One practitioner argues that solver code is becoming a commodity and that verified correctness is what stays scarce.

## 2. Size (order of magnitude only)

| Item | Figure | Source | Confidence |
|---|---|---|---|
| All simulation software | $15.46B (2026) to $28.59B (2031), 13.08% CAGR | Mordor Intelligence | sourced; report mills, treat as rough |
| Optical simulation software | $0.47B (2026) to $1.18B (2035), 11% CAGR | Business Research Insights | sourced; report mill, one number among several that disagree |
| Silicon photonics (devices, not tools) | $2.49B (2025) | InsightAce | sourced; not the tool market |
| Professional photonic design licences | $20k to $100k per licence per year | PhotonDelta blog (ecosystem body) | sourced; not audited |
| COMSOL academic list | about EUR 498 per module per year, EUR 995 perpetual (2018 educational list); about $895 per module on a recent US university price list | COMSOL educational price lists | sourced; commercial prices not retrieved |

Reading: the photonics tool segment is in the hundreds of millions per year, not billions. Academic seats are cheap for the incumbents, so the revenue is in industry licences. **Inferred:** a free solver does not undercut academia on price.

## 3. Who else is there

| Player | What they offer | Relevance to nyx |
|---|---|---|
| Flexcompute | Tidy3D (cloud FDTD, paid), PhotonForge (PIC design, GlobalFoundries and Cadence and Siemens Calibre connectors), FlexAgent MCP, autonomous design loop (reference code on GitHub) | Closest competitor on the agent story. No elastic solver of its own: for acoustics it brings in NGSolve (section 4). |
| COMSOL | FEM incumbent, chatbot in the Model Builder, published guides on connecting coding agents to a COMSOL server | Does the coupled physics today, but a licensed server per session (inferred: limits parallel agent runs). |
| Ansys Lumerical (now under Synopsys) | MODE, FDTD, INTERCONNECT, TCAD | Incumbent for PIC design. |
| Luceda IPKISS, VPIphotonics, Photon Design | Layout and circuit tools | Where the foundry and PDK lock-in lives. |
| femwell (open source, GPL-3.0) | FEM photonic mode solver, thermal solver, inside the gdsfactory ecosystem | The free FEM competitor. Its README lists photonic, electric, thermal and Coulomb solvers, no elastic one (the README says docs lag the code). Last release v0.1.12, Jan 2025; 182 stars. |
| NGSolve (open source) | General FEM library | Used by Flexcompute for the piezoelectric acoustics. A library, not an application: the user or agent writes the model. |
| Meep, MPB, Palace, Elmer | Open-source EM and multiphysics | Need installs and glue code. |
| One lab's own LLM-built solver | A Stanford lab rebuilt a semi-vectorial mode solver in a day with Claude and checked it against a published open solver | Shows the commodity risk for any solver. |

## 4. The gap check (done 2026-11-20)

Optical plus mechanical device design needs optical modes, elastic or piezoelectric modes, and their overlap integrals. What each tool does, from the pages I read:

| Capability | Tidy3D (Flexcompute) | COMSOL | femwell | NGSolve | nyx |
|---|---|---|---|---|---|
| Optical modes | yes, anisotropic mode solver | yes (Wave Optics) | yes | library, you write it | yes (`ems`) |
| Elastic or piezo modes | no; the 17 Nov 2026 article uses NGSolve | yes (Structural Mechanics, Acoustics, MEMS) | none listed in README | yes (used for piezo FEM, 1.1M unknowns) | yes (`sms`: static, eigenfrequency, axisymmetric) |
| Optical-mechanical overlap | agent-written code in a notebook | built into postprocessing (2019 Campinas and Corning story) | optical overlaps only | you write it | by script; the 800 nm SOS intermodal Brillouin case ran this way, awaiting your validation. The `sbs` physics was removed (ADR 24) |
| Agent interface | MCP assistant and an agentic design environment | chatbot, agents via a COMSOL server | none | none | four-tool MCP over a script kernel |
| Layout (GDS) | PhotonForge, foundry connectors | not a layout tool | via gdsfactory | no | none: no mention of GDS, OASIS or KLayout in any source file (grep) |
| Runs where | cloud solver (paid credits) plus local scripts | licensed install | Python, conda | Python or C++ | local, pure Rust; needs the HDF5 library until the static build |

Confidence: the Flexcompute and COMSOL columns are **sourced** from one article each plus product pages. The femwell column is from its README, which says docs lag the code, so a missing solver may exist in the code (unchecked). The nyx column is from the repository and your project notes; I did not re-run anything.

What this means:

- **Flexcompute has shown the whole chain with an agent**, three days before this check. It composes three tools, but a user does not care how many tools sit under the notebook. The remaining weaknesses on their side are two installs, cloud credits and no single validated stack.
- **COMSOL covers it in one tool** with built-in overlap integrals; its barriers are licence cost and a server per session.
- **nyx's edge is the single local stack** (one mesh, one unit system, deterministic, benchmarked against papers) plus no per-run cost, which matters when an agent runs hundreds of sweeps. This is **inferred**, not tested against a user.
- **nyx's weaknesses**: no layout input, no foundry data, Brillouin gain still by script rather than a validated study, and an HDF5 install problem.

Demand:

- **Brillouin devices.** Strong academic activity (a 2025 Twente thesis, a widely cited Brillouin integrated photonics review, thin-film lithium niobate Brillouin lasers), with microwave photonics and low-noise lasers named as targets. The same thesis says there is still a gap between research and deployment. **Inferred:** the commercial design market is small and not yet proven. I found no revenue figure for it.

## 5. What this changes for nyx

1. **Do not sell "AI drives a solver" or "coupled physics" alone.** Both exist. Sell one local, validated stack that an agent can run by the hundred at no per-run cost.
2. **Turn the Brillouin script into a validated study.** The coupling chain (optical mode, elastic mode, overlap, gain) is where Flexcompute leaves the user to write notebook code and COMSOL makes the user click. A study with a published validation case is the product edge. Needs your validation of the 800 nm case first.
3. **Make the benchmark suite a product.** The Tier A cases against papers are the asset. Publish them with the numbers and tolerances. A practitioner's essay argues this is where value moves when solver code becomes easy to replicate (sourced as an opinion, not a measurement).
4. **Add GDS input through gdsfactory.** nyx has none; the open PIC ecosystem builds layouts there and femwell plugs into it. Without it a user must rebuild every layout in nyx's own geometry.
5. **Results format.** The agent loop needs a sweep-level index and per-run provenance, as in `RESULTS_FORMAT_STUDY.md`. Static HDF5 shipping is a precondition for "install and ask Claude".
6. **Be honest about the moat.** Foundry data (PDKs, calibrated models) is the lock-in in this market; nyx has none of it.

## 6. What I would check next

- Look in femwell's code (not only its README) for an elastic solver, and check whether Tidy3D has added one since the 17 Nov article.
- Commercial COMSOL and Lumerical prices (only academic ones found).
- Five or so conversations with people who design Brillouin, optomechanical or MEMS acoustic devices: what they use, what a run costs them, who pays.
- Neighbouring markets that use the same solid-mechanics and mode machinery: MEMS and RF acoustic filters, quantum transduction. Not researched here.

## 7. Sources

- Mordor Intelligence, simulation software market: https://www.mordorintelligence.com/industry-reports/simulation-software-market
- Business Research Insights, optical simulation software market: https://www.businessresearchinsights.com/market-reports/optical-simulation-software-market-115778
- InsightAce, silicon photonics market: https://www.insightaceanalytic.com/report/silicon-photonics-market/3710
- PhotonDelta, photonic chip design software and licence costs: https://www.photondelta.com/blog/what-software-do-you-need-to-design-photonic-chips/
- Flexcompute, FlexAgent MCP (12 Nov 2025) and news list: https://flexcompute.com/news/ai-assisted-simulation-in-tidy3d-ushering-in-a-new-era-of-photonic-design
- Flexcompute, autonomous photonic design: https://github.com/flexcompute/autophotonicdesign and arXiv 2606.00915
- COMSOL, agentic workflows: https://www.comsol.com/blogs/using-agentic-workflows-with-comsol-multiphysics
- Safavi-Naeini, "AI Will Eat Simulation Software": https://quantum-noise.ghost.io/ai-will-eat-simulation-software-and-this-may-be-bad/
- femwell: https://github.com/helgegehring/femwell
- Flexcompute, push-pull acousto-optic modulator on TFLN with Tidy3D, NGSolve and PhotonForge, agent-run (17 Nov 2026): https://flexcompute.com/blog/designing-a-push-pull-acousto-optic-modulator-on-thin-film-lithium-niobate
- COMSOL user story, Brillouin optomechanics at Campinas and Corning (Dec 2019): https://www.comsol.com/story/understanding-brillouin-optomechanical-interaction-with-multiphysics-modeling-80601
- NGSolve: https://ngsolve.org/
- Ye, Brillouin thesis, University of Twente (2025): https://research.utwente.nl/en/publications/stimulated-brillouin-scattering-in-scalable-integrated-photonic-p/
- COMSOL educational price list (university copy): https://techfee.fau.edu/approvedproposals/Download.cfm?sid=514&pid=439
