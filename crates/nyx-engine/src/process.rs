//! Process-wide settings of the program that hosts the engine. Today one: on Windows, opt the process out of power throttling.
//!
//! Windows 11 puts a process it considers background work (a console program, a window that is hidden or not in front, a program started by
//! another one) into "efficiency mode" and runs its threads on the efficiency cores of a hybrid processor. A solve that starts on the
//! performance cores is then moved to the efficiency cores after its first seconds: on an i5-13500 the shifted factorisation of the 3D box
//! cavity took 2.2 s in the first solve of a process and 6.2 s in every later one (`docs/PRODUCTIVITY_STATUS.md`, finding 1). The call below asks
//! Windows to leave this process alone (`PROCESS_POWER_THROTTLING_EXECUTION_SPEED` switched off). It changes which core runs the work and when,
//! never a number: results do not depend on the core or the thread count (ADR 11). Elsewhere it does nothing.

/// Asks Windows not to throttle this process. True when the request was accepted; false when it was refused (an older Windows) or on another
/// system. Call once at the start of `main` of a program that hosts an `Engine`.
pub fn opt_out_of_power_throttling() -> bool {
    #[cfg(windows)]
    {
        #[repr(C)]
        struct PowerThrottlingState {
            version: u32,
            control_mask: u32,
            state_mask: u32,
        }
        #[link(name = "kernel32")]
        extern "system" {
            fn GetCurrentProcess() -> *mut std::ffi::c_void;
            fn SetProcessInformation(
                process: *mut std::ffi::c_void,
                class: i32,
                info: *const PowerThrottlingState,
                size: u32,
            ) -> i32;
        }
        // ProcessPowerThrottling = 4; version 1; control mask 1 = execution speed; state mask 0 = throttling off
        const PROCESS_POWER_THROTTLING: i32 = 4;
        let state = PowerThrottlingState { version: 1, control_mask: 1, state_mask: 0 };
        // SAFETY: the pseudo handle of the current process and a correctly sized, initialised struct that outlives the call.
        unsafe {
            SetProcessInformation(
                GetCurrentProcess(),
                PROCESS_POWER_THROTTLING,
                &state,
                std::mem::size_of::<PowerThrottlingState>() as u32,
            ) != 0
        }
    }
    #[cfg(not(windows))]
    {
        false
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn the_request_is_accepted_on_windows_and_ignored_elsewhere() {
        // a refusal on an old Windows is not a failure of the program; the call must simply not crash
        let accepted = super::opt_out_of_power_throttling();
        if !cfg!(windows) {
            assert!(!accepted);
        }
    }
}
