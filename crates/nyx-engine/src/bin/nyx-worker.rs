//! nyx-worker: solves one frozen project in its own process and writes the results file (productivity plan P4, ADR 86). It is started by the
//! window or the MCP server (`Engine::set_worker`), never by hand: see `nyx_engine::worker` for what it is given and what it answers.
fn main() {
    nyx_engine::process::opt_out_of_power_throttling();
    let args: Vec<String> = std::env::args().skip(1).collect();
    // `nyx-worker queue ...` is `nyx queue ...` for a build that has no window (and for the tests)
    if args.first().map(String::as_str) == Some("queue") {
        std::process::exit(nyx_engine::shared_queue::cli(&args[1..]));
    }
    std::process::exit(nyx_engine::worker::serve(&args));
}
