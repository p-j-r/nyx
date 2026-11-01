Errors are JSON with an "error" code and the data to fix the call.

- unknown_param / bad_value / out_of_range / missing_param / wrong_type: the error lists the valid names, values or range. Fix the call; do not guess.
- unknown_handle: the handle is not in this session (handles live as long as the nyx-mcp process). The error lists the valid handles of that kind. If the session restarted, repeat the steps; the same inputs give the same handles.
- needs_input: not an error. Ask the user the question in "slot"/"why"/"suggest".
- not_guided: the silicon layer is too thin for the wavelength to guide a slab mode. Increase the thickness or shorten the wavelength.
- no_such_mode: the solver returned no mode of that family and order. The error lists the modes found. Try more modes (problem_define n_modes) or a different n_target.
- budget_exceeded: the mesh is larger than this server solves in-process.
- solver_error: femcore refused or failed; the message is its own text.
- A long session may lose the connection to the server. Start a new session and repeat; results are reproducible.
