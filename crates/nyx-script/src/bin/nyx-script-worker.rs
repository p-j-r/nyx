// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The script worker program (productivity plan S4, ADR [89]): `nyx-script-worker subprocess`, started by `nyx-script` through monty-pool. It
//! holds the Monty interpreter and nothing of the engine: every `nyx.<command>(...)` a script makes is sent back to the process that started
//! it as a pause of the script, and that process answers with `Engine::call`. A script that allocates without bound or overflows the stack
//! ends this process, not the engine's.
//!
//! This is `monty-runtime`'s own worker mode (`subprocess.rs`, version 1.6.0, MIT, Copyright Pydantic Services Inc.) kept as a program of this
//! workspace so that building nyx builds its worker and nothing has to be installed beside it. It must be updated with the `monty-*` crates,
//! which are pinned to one version in `Cargo.toml`. The only change: the program has no command line besides the word `subprocess`.

use std::{io, panic, process::ExitCode, thread};

use monty_proto::{
    pb,
    worker::{fatal_error_event, protocol_violation, EventSink, HandleOutcome, SessionBudget, Worker},
    write_frame, FrameError, FrameReader,
};
use monty_types::memory_limit_with_headroom;

/// Bounds the worker's memory and turns an allocation failure into an exit code the pool reports as `MemoryError`. Only a program may
/// install a global allocator, which is why this is a separate program and not a library.
#[global_allocator]
static ALLOC: monty_alloc::LimitedAllocator = monty_alloc::LimitedAllocator;

/// exit code for a protocol violation (the stream with the engine is out of step)
const EX_PROTOCOL: u8 = 76;

/// The interpreter recurses on deep source; a 16 MiB stack leaves room for the nesting limit of the sandbox.
const WORKER_STACK_SIZE: usize = 16 * 1024 * 1024;

fn main() -> ExitCode {
    if std::env::args().nth(1).as_deref() != Some("subprocess") {
        eprintln!("nyx-script-worker is started by nyx, not by hand: `nyx-script-worker subprocess` serves one engine over its standard input and output");
        return ExitCode::from(2);
    }
    thread::Builder::new()
        .name("nyx-script-worker".to_owned())
        .stack_size(WORKER_STACK_SIZE)
        .spawn(serve)
        .expect("failed to spawn the worker thread")
        .join()
        .unwrap_or(ExitCode::from(101))
}

fn serve() -> ExitCode {
    install_panic_hook();
    let mut reader = FrameReader::new(io::stdin().lock());
    let mut worker = Worker::default();
    let mut sink = StdoutSink;
    // the memory baseline is taken on the pristine worker: a session's first request is decoded before its ceiling is armed
    apply_budget(worker.session_budget());

    loop {
        match reader.read::<pb::ClientRequest>() {
            Ok(Some(request)) => {
                let outcome = worker.handle(request, &mut sink);
                apply_budget(worker.session_budget());
                match outcome {
                    Ok(HandleOutcome::Continue) => {}
                    Ok(HandleOutcome::Shutdown) => return ExitCode::SUCCESS,
                    // the worker sent a fatal event (version skew) and cannot go on
                    Ok(HandleOutcome::Fatal) => return ExitCode::from(4),
                    // an oversize event was refused before any byte was written, so the stream is in step and the engine gets a last message
                    Err(FrameError::FrameTooLarge { len, max }) => {
                        fatal(&worker, &mut sink, &format!("response frame of {len} bytes exceeds maximum of {max} bytes"));
                        return ExitCode::from(EX_PROTOCOL);
                    }
                    // writing failed: the engine is gone
                    Err(_) => return ExitCode::from(3),
                }
            }
            // the engine closed our input
            Ok(None) => return ExitCode::SUCCESS,
            // the frame arrived whole but its payload did not decode: answer with an error that ends the turn and go on
            Err(FrameError::Decode(err)) => {
                if sink.send(&protocol_violation(&format!("malformed request: {err}"))).is_err() {
                    return ExitCode::from(3);
                }
            }
            // the stream is out of step: nothing can be recovered
            Err(err) => {
                fatal(&worker, &mut sink, &format!("malformed request frame: {err}"));
                return ExitCode::from(EX_PROTOCOL);
            }
        }
    }
}

fn apply_budget(budget: SessionBudget) {
    let hard_memory_limit = memory_limit_with_headroom(budget.max_memory, budget.type_check);
    monty_alloc::set_hard_limit(hard_memory_limit).expect("the allocator of this program must be monty-alloc's");
}

struct StdoutSink;

impl EventSink for StdoutSink {
    fn send(&mut self, event: &pb::WorkerEvent) -> Result<(), FrameError> {
        write_frame(&mut io::stdout(), event)
    }

    fn session_opened(&mut self, budget: SessionBudget) {
        apply_budget(budget);
    }
}

fn fatal(worker: &Worker, sink: &mut impl EventSink, message: &str) {
    eprintln!("nyx-script-worker fatal error: {message}");
    let _ = sink.send(&worker.fatal_event(message));
}

/// A panic is reported to the engine as a fatal event first (the engine then knows it was a crash and says so), then handled as usual.
fn install_panic_hook() {
    let default_hook = panic::take_hook();
    panic::set_hook(Box::new(move |info| {
        let _ = write_frame(&mut io::stdout(), &fatal_error_event(&format!("worker panicked: {info}")));
        default_hook(info);
    }));
}
