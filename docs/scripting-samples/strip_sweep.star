# strip_sweep.star -- SAMPLE ONLY: nyx has no Starlark binding yet. Written to compare the language with strip_sweep.rhai.
#
# Assumed binding (what the host would register; generated from Engine::catalog()):
#   nyx.call(name, **args)          any engine command, the same JSON as an MCP tool call
#   nyx.add_rect(name=..., min=..., max=...), nyx.set_material(...), ...   one function per catalog command, keyword arguments
#   nyx.try_call(name, **args)      -> struct(ok, value, error)   Starlark has no try/except, so recoverable errors need this
# An engine error raised by nyx.<command> stops the whole script (like fail()).
# Starlark itself has no file, network, clock or random access; the host gives none either.
#   fmt(x, digits)                  fixed-point text of a number: Starlark's `%` operator has no precision ("%.5f" is an error), so the host adds this (the spike found it)
#
# Workflow: a silicon strip on oxide at 1550 nm.
#   1. build the model with a parameter `w` (strip width)
#   2. mesh convergence on TE0: refine until neff changes by less than 2e-3 (1e-4 needs a 10 nm mesh, 40 minutes: measured by the spike)
#   3. sweep the width at the converged mesh size, return a small table

C0 = 299792458.0
LAMBDA = 1.55e-6
TOL = 2e-3

def build():
    nyx.project_new(dim = 2, title = "SOI strip width sweep")
    nyx.set_param(name = "w", value = 500e-9)
    nyx.set_param(name = "t", value = 220e-9)

    nyx.add_rect(name = "box", min = [-2.0e-6, -2.0e-6], max = [2.0e-6, 0.0])
    nyx.add_rect(name = "clad", min = [-2.0e-6, 0.0], max = [2.0e-6, 2.0e-6])
    nyx.add_rect(name = "core", min = ["-w/2", 0.0], max = ["w/2", "t"])  # fields may be parameter expressions
    nyx.add_overlay(regions = [
        {"name": "BOX", "input": "box"},
        {"name": "Clad", "input": "clad"},
        {"name": "Si", "input": "core"},  # painter order: the core wins
    ])

    nyx.add_physics(
        physics = "ems",
        study = "modes",
        study_params = {"freq": C0 / LAMBDA, "n_target": 2.9, "n_modes": 4},
    )
    for region, material in [("BOX", "silica"), ("Clad", "silica"), ("Si", "silicon")]:
        nyx.set_material(physics = "ems", region = region, material = material)

def solve_modes(solve = True):
    if solve:
        nyx.solve_run(physics = "ems", wait = True)
    res = nyx.result_info(physics = "ems")[0]
    if res["stale"]:
        fail("result is stale")
    pol = nyx.analyse(analysis = "polarization", physics = "ems")

    # shape of the polarization rows: { label, value:[re,im] (the effective index), report:{ class: "quasi-TE (x)", te_fraction, ... } }
    te = None
    tm = None
    for row in pol["rows"]:
        cls = row["report"]["class"]
        if cls.startswith("quasi-TE") and te == None:
            te = row["value"][0]
        if cls.startswith("quasi-TM") and tm == None:
            tm = row["value"][0]
    return {"te0": te, "tm0": tm, "dofs": res["dofs"], "seconds": res["seconds"]}

def converge(h0, h_min):
    # Standard Starlark has no `while` and no recursion: every loop runs over a finite list,
    # so the script provably ends. The refinement ladder is written out instead.
    ladder = []
    h = h0
    for _ in range(8):
        if h < h_min:
            break
        ladder.append(h)
        h = h / 2

    prev = None
    for i, h in enumerate(ladder):
        nyx.set_mesh(h = h, adapt_scale = 0.4)
        cur = solve_modes()
        if prev != None:
            d = abs(cur["te0"] - prev["te0"])
            print("h = " + fmt(h * 1e9, 1) + " nm  TE0 = " + fmt(cur["te0"], 6) + "  change = " + fmt(d, 6) + "  dofs = " + str(cur["dofs"]))
            if d < TOL:
                return ladder[i - 1]  # the coarser of the last pair already met the tolerance
        prev = cur
    fail("TE0 not converged to " + fmt(TOL, 6) + " above h = " + fmt(h_min * 1e9, 0) + " nm")

def sweep(h):
    nyx.set_mesh(h = h, adapt_scale = 0.4)
    rows = []
    for w_nm in range(400, 701, 100):
        nyx.set_param(name = "w", value = w_nm * 1e-9)
        r = nyx.try_call("solve_run", physics = "ems", wait = True)  # keep the other points if one fails
        if not r.ok:
            rows.append({"w_nm": w_nm, "error": r.error})
            continue
        m = solve_modes(solve = False)  # already solved by try_call
        rows.append({"w_nm": w_nm, "te0": m["te0"], "tm0": m["tm0"]})
    return rows

def main():
    build()
    h = converge(200e-9, 25e-9)
    rows = sweep(h)

    # Only what is printed or returned goes back to the agent; fields stay in the engine.
    print("converged h = " + fmt(h * 1e9, 1) + " nm")
    print("w_nm  TE0       TM0")
    for r in rows:
        if "error" in r:
            print(str(r["w_nm"]) + "   failed: " + r["error"])
        else:
            print(str(r["w_nm"]) + "   " + fmt(r["te0"], 5) + "   " + fmt(r["tm0"], 5))
    return rows

result = main()  # the host returns `result` to the caller as JSON
