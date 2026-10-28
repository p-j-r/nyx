//! XML for the model description: any serde value (the whole project: geometry, materials, boundary
//! conditions, studies) is written as readable XML with `quick-xml` and read back exactly.
//!
//! Encoding of a value, always with a `t` attribute naming its type:
//!   object  `<Name t="o">` with one child per field, named after the key (or `<entry key="1">`
//!           when the key is not a valid XML name, e.g. a domain number)
//!   array   `<Name t="a">` with `<item>` children
//!   string  `<Name t="s">text</Name>`   number `t="n"`   boolean `t="b"`   null `t="z"`

use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::{Reader, Writer};
use serde_json::{Map, Number, Value};

fn valid_name(k: &str) -> bool {
    let mut c = k.chars();
    matches!(c.next(), Some(ch) if ch.is_ascii_alphabetic() || ch == '_') && c.all(|ch| ch.is_ascii_alphanumeric() || ch == '_' || ch == '-')
        && !k.to_ascii_lowercase().starts_with("xml")
        && k != "entry"
        && k != "item"
}

fn write_value(w: &mut Writer<Vec<u8>>, name: &str, key_attr: Option<&str>, v: &Value) -> Result<(), String> {
    let mut start = BytesStart::new(name);
    if let Some(k) = key_attr {
        start.push_attribute(("key", k));
    }
    let t = match v {
        Value::Object(_) => "o",
        Value::Array(_) => "a",
        Value::String(_) => "s",
        Value::Number(_) => "n",
        Value::Bool(_) => "b",
        Value::Null => "z",
    };
    start.push_attribute(("t", t));
    let err = |e: std::io::Error| e.to_string();
    w.write_event(Event::Start(start)).map_err(err)?;
    match v {
        Value::Object(m) => {
            for (k, x) in m {
                if valid_name(k) {
                    write_value(w, k, None, x)?;
                } else {
                    write_value(w, "entry", Some(k), x)?;
                }
            }
        }
        Value::Array(a) => {
            for x in a {
                write_value(w, "item", None, x)?;
            }
        }
        Value::String(s) => w.write_event(Event::Text(BytesText::new(s))).map_err(err)?,
        Value::Number(n) => w.write_event(Event::Text(BytesText::new(&n.to_string()))).map_err(err)?,
        Value::Bool(b) => w.write_event(Event::Text(BytesText::new(if *b { "true" } else { "false" }))).map_err(err)?,
        Value::Null => {}
    }
    w.write_event(Event::End(BytesEnd::new(name))).map_err(err)
}

/// A JSON value as indented XML under the root element `root`.
pub fn value_to_xml(root: &str, v: &Value) -> Result<String, String> {
    let mut w = Writer::new_with_indent(Vec::new(), b' ', 2);
    w.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None))).map_err(|e| e.to_string())?;
    write_value(&mut w, root, None, v)?;
    String::from_utf8(w.into_inner()).map_err(|e| e.to_string())
}

struct Frame {
    name: String,
    key: Option<String>,
    t: String,
    text: String,
    obj: Map<String, Value>,
    arr: Vec<Value>,
}

fn open_frame(e: &BytesStart) -> Result<Frame, String> {
    let name = String::from_utf8_lossy(e.name().as_ref()).into_owned();
    let (mut t, mut key) = (String::new(), None);
    for a in e.attributes().flatten() {
        let v = a.normalized_value(quick_xml::XmlVersion::Explicit1_0).map_err(|x| x.to_string())?.into_owned();
        match a.key.as_ref() {
            b"t" => t = v,
            b"key" => key = Some(v),
            _ => {}
        }
    }
    if t.is_empty() {
        return Err(format!("element <{name}> has no type attribute: not a nyx model file"));
    }
    Ok(Frame { name, key, t, text: String::new(), obj: Map::new(), arr: vec![] })
}

fn finish(f: &Frame) -> Result<Value, String> {
    Ok(match f.t.as_str() {
        "o" => Value::Object(f.obj.clone()),
        "a" => Value::Array(f.arr.clone()),
        "s" => Value::String(f.text.clone()),
        "n" => {
            let s = f.text.trim();
            if let Ok(i) = s.parse::<i64>() {
                Value::Number(i.into())
            } else if let Ok(u) = s.parse::<u64>() {
                Value::Number(u.into())
            } else {
                let x: f64 = s.parse().map_err(|_| format!("bad number '{s}' in <{}>", f.name))?;
                Value::Number(Number::from_f64(x).ok_or("non-finite number")?)
            }
        }
        "b" => Value::Bool(f.text.trim() == "true"),
        "z" => Value::Null,
        other => return Err(format!("unknown type '{other}' in <{}>", f.name)),
    })
}

fn attach(stack: &mut [Frame], result: &mut Option<(String, Value)>, f: Frame, val: Value) -> Result<(), String> {
    match stack.last_mut() {
        Some(p) => match p.t.as_str() {
            "o" => {
                let k = if f.name == "entry" { f.key.ok_or("entry without key")? } else { f.name };
                p.obj.insert(k, val);
                Ok(())
            }
            "a" => {
                p.arr.push(val);
                Ok(())
            }
            _ => Err(format!("<{}> is not a container", p.name)),
        },
        None => {
            *result = Some((f.name, val));
            Ok(())
        }
    }
}

/// Reads XML written by `value_to_xml`; returns the root element name and its value.
pub fn xml_to_value(xml: &str) -> Result<(String, Value), String> {
    let mut reader = Reader::from_str(xml);
    let mut stack: Vec<Frame> = vec![];
    let mut result: Option<(String, Value)> = None;
    loop {
        match reader.read_event().map_err(|e| format!("XML error at byte {}: {e}", reader.buffer_position()))? {
            Event::Start(e) => stack.push(open_frame(&e)?),
            Event::Empty(e) => {
                let f = open_frame(&e)?;
                let val = finish(&f)?;
                attach(&mut stack, &mut result, f, val)?;
            }
            Event::Text(t) => {
                if let Some(f) = stack.last_mut() {
                    f.text.push_str(&t.decode().map_err(|x| x.to_string())?);
                }
            }
            Event::GeneralRef(r) => {
                if let Some(f) = stack.last_mut() {
                    let name = String::from_utf8_lossy(&r).into_owned();
                    f.text.push_str(match name.as_str() {
                        "amp" => "&",
                        "lt" => "<",
                        "gt" => ">",
                        "quot" => "\"",
                        "apos" => "'",
                        other => return Err(format!("unknown entity &{other};")),
                    });
                }
            }
            Event::End(_) => {
                let f = stack.pop().ok_or("unbalanced XML")?;
                let val = finish(&f)?;
                attach(&mut stack, &mut result, f, val)?;
            }
            Event::Eof => break,
            _ => {}
        }
    }
    result.ok_or_else(|| "empty XML document".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_awkward_values() {
        let v = serde_json::json!({
            "name": "a < b & \"c\"", "empty": "", "n": 1.5e-9, "i": -3, "big": 1e300, "flag": true, "nothing": null,
            "list": [1, [2, 3], {"x": []}, "s"],
            "1": {"E": 2e11}, "weird key": 7, "item": 1, "entry": 2, "xml1": 3
        });
        let x = value_to_xml("model", &v).unwrap();
        let (root, back) = xml_to_value(&x).unwrap();
        assert_eq!(root, "model");
        assert_eq!(back, v, "{x}");
    }

    #[test]
    fn rejects_foreign_xml() {
        assert!(xml_to_value("<a><b/></a>").is_err());
        assert!(xml_to_value("not xml at all <").is_err() || !xml_to_value("not xml at all <").unwrap_err().is_empty());
    }
}
