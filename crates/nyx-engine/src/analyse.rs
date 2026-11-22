// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Analyses of results as engine commands: a result is solved once, an analysis is a function of it and of the mesh it was solved on
//! (`nyx-analysis`, ADR 56). The window draws the typed `ModeAnalysis`; the MCP server gets the same through `analyse`, as plain JSON.

use crate::engine::Engine;
use nyx_analysis::{
    Analysis, CapacitanceMatrixAnalysis, EprAnalysis, HamiltonianAnalysis, InductanceMatrixAnalysis, ModeAnalysis,
    NetworkAnalysis, PolarizationAnalysis, TlsParticipationAnalysis,
};
use nyx_physics::Error;
use serde_json::{json, Value};

/// The analyses of this build.
pub fn analyses() -> Vec<Box<dyn Analysis>> {
    vec![
        Box::new(PolarizationAnalysis),
        Box::new(CapacitanceMatrixAnalysis),
        Box::new(InductanceMatrixAnalysis),
        Box::new(EprAnalysis),
        Box::new(HamiltonianAnalysis),
        Box::new(TlsParticipationAnalysis),
        Box::new(NetworkAnalysis),
    ]
}

impl Engine {
    /// The mesh a result was solved on, or why there is none (the geometry or mesh settings changed since).
    pub(crate) fn analysis_inputs(
        &self,
        physics: &str,
    ) -> Result<(std::sync::Arc<crate::meshing::MeshedModel>, std::sync::Arc<nyx_physics::SolveResult>), Error> {
        let stored = self
            .results
            .get(physics)
            .ok_or_else(|| Error::invalid(format!("there is no result for '{physics}' yet; solve first")))?;
        let model = self
            .model
            .as_ref()
            .map(|(_, m)| m.clone())
            .ok_or_else(|| Error::invalid("there is no mesh; solve first"))?;
        Ok((model, stored.result.clone()))
    }

    /// The polarization analysis of the modes of `physics`, typed (what the Analysis tab draws).
    pub fn polarization(&self, physics: &str) -> Result<ModeAnalysis, Error> {
        let (model, result) = self.analysis_inputs(physics)?;
        nyx_analysis::polarization(model.mesh(), &result)
    }

    /// Whether the polarization applies to the latest result of `physics`.
    pub fn polarization_applies(&self, physics: &str) -> bool {
        self.results.get(physics).is_some_and(|s| nyx_analysis::polarization_applies(&s.result))
    }

    /// Lists the analyses and says which apply to which result.
    pub fn list_analyses(&self) -> Value {
        json!(analyses()
            .iter()
            .map(|a| json!({
                "analysis": a.id(),
                "name": a.name(),
                "applies_to": self.results.iter().filter(|(_, s)| a.applies_to(&s.result)).map(|(k, _)| k.clone()).collect::<Vec<_>>(),
                "stale": self.results.iter().filter(|(_, s)| a.applies_to(&s.result)).map(|(k, _)| (k.clone(), self.is_stale(k))).collect::<std::collections::BTreeMap<_, _>>(),
            }))
            .collect::<Vec<_>>())
    }

    /// Runs one analysis on the latest result of `physics` and returns its JSON.
    pub fn analyse(&self, analysis: &str, physics: &str) -> Result<Value, Error> {
        self.analyse_with(analysis, physics, &Value::Null)
    }

    /// Runs one analysis on the latest result of `physics` with optional `params` and returns its JSON.
    pub fn analyse_with(&self, analysis: &str, physics: &str, params: &Value) -> Result<Value, Error> {
        let all = analyses();
        let a = all.iter().find(|a| a.id() == analysis).ok_or_else(|| {
            Error::invalid(format!(
                "there is no analysis '{analysis}' (they are: {})",
                all.iter().map(|a| a.id()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        let (model, result) = self.analysis_inputs(physics)?;
        if !a.applies_to(&result) {
            return Err(Error::invalid(format!(
                "the {analysis} analysis does not apply to the {} result of '{physics}'",
                result.study
            )));
        }
        let mut v = a.run_with_params(model.mesh(), &result, params)?;
        if let Some(o) = v.as_object_mut() {
            o.insert("analysis".into(), json!(analysis));
            o.insert("physics".into(), json!(physics));
            o.insert("stale".into(), json!(self.is_stale(physics)));
            o.insert("age_s".into(), json!(self.result_age(physics)));
            o.insert("validated".into(), self.validated_json(physics));
        }
        Ok(v)
    }
}
