//! Adaptation stack. See the milestone table in `docs/STATUS.md`.

pub mod drive;
pub mod drive3;
pub mod electro;
pub mod field;
pub mod field3;
pub mod geom;
pub mod geom3;
pub mod maxwell;
pub mod metric;
pub mod ops;
pub mod ops3;
pub mod topo;
pub mod topo3;
pub mod estimate;
pub mod estimate3;
pub mod timer;
pub mod work;

/// What the estimators read of a result: the plot fields of every solution as `(flat name, values per dof)`. Implemented for the result
/// of a physics crate in `nyx-physics`.
pub trait PlotFields {
    fn plot_fields(&self) -> Vec<Vec<(String, Vec<f64>)>>;
}
