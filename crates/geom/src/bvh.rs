// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Zero-dependency flat `Vec`-backed 2D AABB bounding volume hierarchy (`AabbBvh2`) and complex
//! spatial index (`ComplexBvh`) for $O(\log N)$ region, loop, and edge queries.

use crate::math::{Aabb, V2};

const LEAF_CAP: usize = 8;

/// 2D axis-aligned bounding box in the planar curve `(u, v)` frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bbox2 {
    pub lo: V2,
    pub hi: V2,
}

impl Bbox2 {
    pub const EMPTY: Bbox2 = Bbox2 { lo: [f64::INFINITY; 2], hi: [f64::NEG_INFINITY; 2] };

    #[inline]
    pub fn new(lo: V2, hi: V2) -> Self {
        Self { lo: [lo[0].min(hi[0]), lo[1].min(hi[1])], hi: [lo[0].max(hi[0]), lo[1].max(hi[1])] }
    }

    #[inline]
    pub fn from_aabb_xy(b: &Aabb) -> Self {
        if b.is_empty() {
            Self::EMPTY
        } else {
            Self { lo: [b.lo[0], b.lo[1]], hi: [b.hi[0], b.hi[1]] }
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.lo[0] > self.hi[0] || self.lo[1] > self.hi[1]
    }

    #[inline]
    pub fn grow(&mut self, p: V2) {
        self.lo[0] = self.lo[0].min(p[0]);
        self.lo[1] = self.lo[1].min(p[1]);
        self.hi[0] = self.hi[0].max(p[0]);
        self.hi[1] = self.hi[1].max(p[1]);
    }

    #[inline]
    pub fn union(&self, o: &Bbox2) -> Bbox2 {
        if self.is_empty() {
            *o
        } else if o.is_empty() {
            *self
        } else {
            Bbox2 {
                lo: [self.lo[0].min(o.lo[0]), self.lo[1].min(o.lo[1])],
                hi: [self.hi[0].max(o.hi[0]), self.hi[1].max(o.hi[1])],
            }
        }
    }

    #[inline]
    pub fn pad(&self, d: f64) -> Bbox2 {
        if self.is_empty() {
            *self
        } else {
            Bbox2 { lo: [self.lo[0] - d, self.lo[1] - d], hi: [self.hi[0] + d, self.hi[1] + d] }
        }
    }

    #[inline]
    pub fn contains_point(&self, p: V2, pad: f64) -> bool {
        p[0] >= self.lo[0] - pad && p[0] <= self.hi[0] + pad && p[1] >= self.lo[1] - pad && p[1] <= self.hi[1] + pad
    }

    #[inline]
    pub fn overlaps(&self, o: &Bbox2, pad: f64) -> bool {
        self.lo[0] <= o.hi[0] + pad
            && o.lo[0] <= self.hi[0] + pad
            && self.lo[1] <= o.hi[1] + pad
            && o.lo[1] <= self.hi[1] + pad
    }

    /// True if this box overlaps the horizontal ray `[p[0], +inf) x {p[1]}` within `pad`.
    #[inline]
    pub fn overlaps_ray_x(&self, p: V2, pad: f64) -> bool {
        self.hi[0] + pad >= p[0] && self.lo[1] - pad <= p[1] && self.hi[1] + pad >= p[1]
    }

    #[inline]
    pub fn center(&self) -> V2 {
        [0.5 * (self.lo[0] + self.hi[0]), 0.5 * (self.lo[1] + self.hi[1])]
    }
}

/// One node in a flat `Vec`-backed 2D AABB BVH.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BvhNode2 {
    pub bbox: Bbox2,
    /// For internal nodes (`count == 0`): left child index in `nodes`.
    /// For leaf nodes (`count > 0`): start offset in `items`.
    pub left_or_start: u32,
    /// For internal nodes (`count == 0`): right child index in `nodes`.
    pub right: u32,
    /// Number of items in a leaf (`0` for internal nodes).
    pub count: u32,
}

/// Flat, deterministic 2D AABB bounding volume hierarchy.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AabbBvh2 {
    pub nodes: Vec<BvhNode2>,
    pub items: Vec<u32>,
    pub item_boxes: Vec<Bbox2>,
}

impl AabbBvh2 {
    /// Builds a deterministic top-down median-split 2D BVH over `item_boxes`.
    pub fn build(item_boxes: Vec<Bbox2>) -> Self {
        let mut work: Vec<(u32, Bbox2, V2)> = item_boxes
            .iter()
            .enumerate()
            .filter(|(_, b)| !b.is_empty())
            .map(|(i, &b)| (i as u32, b, b.center()))
            .collect();
        if work.is_empty() {
            return Self { nodes: Vec::new(), items: Vec::new(), item_boxes };
        }
        let mut nodes = Vec::with_capacity(2 * (work.len() / LEAF_CAP + 1));
        let mut items = Vec::with_capacity(work.len());
        Self::build_rec(&mut work, &mut nodes, &mut items);
        Self { nodes, items, item_boxes }
    }

    fn build_rec(work: &mut [(u32, Bbox2, V2)], nodes: &mut Vec<BvhNode2>, items: &mut Vec<u32>) -> u32 {
        let node_idx = nodes.len() as u32;
        let mut bbox = Bbox2::EMPTY;
        let mut cbox = Bbox2::EMPTY;
        for &(_, ref b, c) in work.iter() {
            bbox = bbox.union(b);
            cbox.grow(c);
        }
        if work.len() <= LEAF_CAP {
            work.sort_unstable_by_key(|x| x.0);
            let start = items.len() as u32;
            let count = work.len() as u32;
            for &(id, _, _) in work.iter() {
                items.push(id);
            }
            nodes.push(BvhNode2 { bbox, left_or_start: start, right: 0, count });
            return node_idx;
        }
        nodes.push(BvhNode2 { bbox, left_or_start: 0, right: 0, count: 0 });
        let axis = if (cbox.hi[0] - cbox.lo[0]) >= (cbox.hi[1] - cbox.lo[1]) { 0 } else { 1 };
        let mid = work.len() / 2;
        work.select_nth_unstable_by(mid, |a, b| a.2[axis].total_cmp(&b.2[axis]).then(a.0.cmp(&b.0)));
        let (left_slice, right_slice) = work.split_at_mut(mid);
        let left_child = Self::build_rec(left_slice, nodes, items);
        let right_child = Self::build_rec(right_slice, nodes, items);
        nodes[node_idx as usize].left_or_start = left_child;
        nodes[node_idx as usize].right = right_child;
        node_idx
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    #[inline]
    fn traverse(&self, mut matches: impl FnMut(&Bbox2) -> bool, mut visit: impl FnMut(u32)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut inline_stack = [0u32; 64];
        let mut sp = 1usize;
        let mut spill: Vec<u32> = Vec::new();
        let mut spilled = false;
        loop {
            let idx = if spilled {
                let Some(v) = spill.pop() else { break };
                v
            } else {
                if sp == 0 {
                    break;
                }
                sp -= 1;
                inline_stack[sp]
            };
            let node = &self.nodes[idx as usize];
            if !matches(&node.bbox) {
                continue;
            }
            if node.count > 0 {
                let start = node.left_or_start as usize;
                let end = start + node.count as usize;
                for &id in &self.items[start..end] {
                    if matches(&self.item_boxes[id as usize]) {
                        visit(id);
                    }
                }
            } else if spilled {
                spill.push(node.right);
                spill.push(node.left_or_start);
            } else if sp + 2 <= inline_stack.len() {
                inline_stack[sp] = node.right;
                inline_stack[sp + 1] = node.left_or_start;
                sp += 2;
            } else {
                spill.reserve(sp + 32);
                spill.extend_from_slice(&inline_stack[..sp]);
                spill.push(node.right);
                spill.push(node.left_or_start);
                spilled = true;
            }
        }
    }

    /// Visits all item indices whose 2D AABB (padded by `pad`) contains `p`.
    #[inline]
    pub fn query_point(&self, p: V2, pad: f64, visit: impl FnMut(u32)) {
        self.traverse(|bb| bb.contains_point(p, pad), visit);
    }

    /// Visits all item indices whose 2D AABB overlaps `target` within `pad`.
    #[inline]
    pub fn query_box(&self, target: &Bbox2, pad: f64, visit: impl FnMut(u32)) {
        self.traverse(|bb| bb.overlaps(target, pad), visit);
    }

    /// Visits all item indices whose 2D AABB overlaps the horizontal ray `[p[0], +inf) x {p[1]}` within `pad`.
    #[inline]
    pub fn query_ray_x(&self, p: V2, pad: f64, visit: impl FnMut(u32)) {
        self.traverse(|bb| bb.overlaps_ray_x(p, pad), visit);
    }
}

/// Spatial acceleration hierarchy attached to a 2D [`crate::complex::Complex`] (in the 2D curve `(u, v)` frame).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComplexBvh {
    pub regions: AabbBvh2,
    pub loops: AabbBvh2,
    pub edges: AabbBvh2,
    /// Per-loop BVH over `0..loop.edges.len()` for loops with `>= 32` edges.
    pub loop_edge_bvhs: Vec<Option<AabbBvh2>>,
}
