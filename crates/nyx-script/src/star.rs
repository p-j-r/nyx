//! The Starlark binding (S4 spike): `nyx.<command>(**args)` for every command of the catalog, `nyx.call(name, **args)` and `nyx.try_call(name, **args)`.
//!
//! The functions are a generated prelude (Starlark source) over two native functions, so the binding follows `Engine::catalog()` with no hand-written
//! list. An engine error stops the script like `fail()`; `try_call` returns `struct(ok, value, error)` for the recoverable case.

use std::cell::{Cell, RefCell};

use nyx_engine::Engine;
use serde_json::Value;
use starlark::any::ProvidesStaticType;
use starlark::environment::{FrozenModule, Globals, GlobalsBuilder, LibraryExtension, Module};
use starlark::eval::Evaluator;
use starlark::starlark_module;
use starlark::PrintHandler;
use starlark::syntax::{AstModule, Dialect};
use starlark::values::list::AllocList;
use starlark::values::dict::AllocDict;
use starlark::values::structs::AllocStruct;
use starlark::values::float::UnpackFloat;
use starlark::values::{Heap, Value as SV};

use crate::{Limits, Outcome};

#[derive(ProvidesStaticType)]
struct Host<'a> {
    engine: RefCell<&'a mut Engine>,
    calls: Cell<u64>,
}

#[derive(Default)]
struct Printed(RefCell<Vec<String>>);

impl PrintHandler for Printed {
    fn println(&self, text: &str) -> starlark::Result<()> {
        self.0.borrow_mut().push(text.to_string());
        Ok(())
    }
}

fn to_star<'v>(j: &Value, heap: Heap<'v>) -> SV<'v> {
    match j {
        Value::Null => SV::new_none(),
        Value::Bool(b) => SV::new_bool(*b),
        Value::Number(n) => match n.as_i64() {
            Some(i) => heap.alloc(i),
            None => heap.alloc(n.as_f64().unwrap_or(f64::NAN)),
        },
        Value::String(s) => heap.alloc(s.as_str()),
        Value::Array(a) => heap.alloc(AllocList(a.iter().map(|x| to_star(x, heap)))),
        Value::Object(o) => heap.alloc(AllocDict(o.iter().map(|(k, v)| (k.as_str(), to_star(v, heap))))),
    }
}

fn engine_call<'v>(host: &Host, name: &str, args: SV<'v>) -> Result<Value, String> {
    let args = args.to_json_value().map_err(|e| format!("{name}: arguments are not JSON: {e}"))?;
    host.calls.set(host.calls.get() + 1);
    host.engine.borrow_mut().call(name, &args).map_err(|e| e.to_string())
}

fn host<'v, 'a, 'e>(eval: &Evaluator<'v, 'a, 'e>) -> anyhow::Result<&'a Host<'e>> {
    let any = eval.extra.ok_or_else(|| anyhow::anyhow!("no engine is attached to this script"))?;
    any.downcast_ref::<Host>().ok_or_else(|| anyhow::anyhow!("no engine is attached to this script"))
}

#[starlark_module]
fn natives(b: &mut GlobalsBuilder) {
    fn fmt(x: UnpackFloat, digits: i32) -> anyhow::Result<String> {
        Ok(crate::fmt_fixed(x.0, digits as i64))
    }
    fn _nyx_call<'v>(name: &str, args: SV<'v>, eval: &mut Evaluator<'v, '_, '_>) -> anyhow::Result<SV<'v>> {
        let h = host(eval)?;
        let v = engine_call(h, name, args).map_err(|e| anyhow::anyhow!(e))?;
        Ok(to_star(&v, eval.heap()))
    }
    fn _nyx_try_call<'v>(name: &str, args: SV<'v>, eval: &mut Evaluator<'v, '_, '_>) -> anyhow::Result<SV<'v>> {
        let h = host(eval)?;
        let heap = eval.heap();
        Ok(match engine_call(h, name, args) {
            Ok(v) => heap.alloc(AllocStruct([("ok", SV::new_bool(true)), ("value", to_star(&v, heap)), ("error", SV::new_none())])),
            Err(e) => heap.alloc(AllocStruct([("ok", SV::new_bool(false)), ("value", SV::new_none()), ("error", heap.alloc(e.as_str()))])),
        })
    }
}

/// The Starlark source of the `nyx` struct, one function per command.
pub fn prelude(engine: &Engine) -> String {
    let mut s = String::from("def _mk(n):\n    def f(**kw):\n        return _nyx_call(n, kw)\n    return f\n\n");
    s.push_str("def _call(name, **kw):\n    return _nyx_call(name, kw)\n\ndef _try_call(name, **kw):\n    return _nyx_try_call(name, kw)\n\n");
    s.push_str("nyx = struct(\n    call = _call,\n    try_call = _try_call,\n");
    for c in engine.catalog() {
        s.push_str(&format!("    {} = _mk(\"{}\"),\n", c.name, c.name));
    }
    s.push_str(")\n");
    s
}

fn globals() -> Globals {
    GlobalsBuilder::extended_by(&[LibraryExtension::StructType, LibraryExtension::Print]).with(natives).build()
}

fn load_prelude(src: &str) -> Result<FrozenModule, String> {
    let ast = AstModule::parse("nyx_prelude", src.to_string(), &Dialect::Standard).map_err(|e| e.to_string())?;
    Module::with_temp_heap(|module| {
        {
            let mut eval = Evaluator::new(&module);
            eval.eval_module(ast, &globals()).map_err(|e| e.to_string())?;
        }
        module.freeze().map_err(|e| format!("{e:?}"))
    })
}

/// Parses the script and nothing else (the static check of stage 1 reuses this).
pub fn parse(name: &str, src: &str) -> Result<(), String> {
    AstModule::parse(name, src.to_string(), &Dialect::Standard).map(|_| ()).map_err(|e| e.to_string())
}

/// Runs a script. Its value is the global `result` when it defines one.
pub fn run(engine: &mut Engine, name: &str, src: &str, limits: &Limits) -> Result<Outcome, String> {
    let pre = load_prelude(&prelude(engine))?;
    let ast = AstModule::parse(name, src.to_string(), &Dialect::Standard).map_err(|e| e.to_string())?;
    let printed = Printed::default();
    let host = Host { engine: RefCell::new(engine), calls: Cell::new(0) };
    let value = Module::with_temp_heap(|module| -> Result<Value, String> {
        module.import_public_symbols(&pre);
        let mut eval = Evaluator::new(&module);
        eval.extra = Some(&host);
        eval.set_print_handler(&printed);
        eval.set_max_tick_count(limits.steps).map_err(|e| e.to_string())?;
        eval.set_max_heap_size(limits.heap_bytes).map_err(|e| e.to_string())?;
        eval.set_max_callstack_size(limits.depth).map_err(|e| e.to_string())?;
        eval.eval_module(ast, &globals()).map_err(|e| e.to_string())?;
        Ok(match module.get("result") {
            Some(r) => r.to_json_value().map_err(|e| e.to_string())?,
            None => Value::Null,
        })
    })?;
    Ok(Outcome { value, printed: printed.0.into_inner(), calls: host.calls.get() })
}
