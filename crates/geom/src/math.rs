//! Small vector, box and similarity-transform helpers (2D and 3D, `f64`, plain arrays).

pub type V2 = [f64; 2];
pub type V3 = [f64; 3];

pub const TAU: f64 = std::f64::consts::TAU;

#[inline]
pub fn add2(a: V2, b: V2) -> V2 {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline]
pub fn sub2(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
pub fn mul2(a: V2, s: f64) -> V2 {
    [a[0] * s, a[1] * s]
}
#[inline]
pub fn dot2(a: V2, b: V2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
#[inline]
pub fn cross2(a: V2, b: V2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
#[inline]
pub fn norm2(a: V2) -> f64 {
    a[0].hypot(a[1])
}
#[inline]
pub fn dist2(a: V2, b: V2) -> f64 {
    norm2(sub2(a, b))
}
/// Rotation by +90 degrees (the left normal of a direction).
#[inline]
pub fn perp(a: V2) -> V2 {
    [-a[1], a[0]]
}

#[inline]
pub fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub fn add3(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub fn mul3(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub fn dot3(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub fn cross3(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
pub fn norm3(a: V3) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
pub fn unit3(a: V3) -> V3 {
    let n = norm3(a);
    if n > 0.0 {
        mul3(a, 1.0 / n)
    } else {
        a
    }
}
#[inline]
pub fn to3(p: V2) -> V3 {
    [p[0], p[1], 0.0]
}

/// Angle of `a` relative to `base`, wrapped into `[base, base + TAU)`.
#[inline]
pub fn wrap_from(t: f64, base: f64) -> f64 {
    base + (t - base).rem_euclid(TAU)
}

/// Axis-aligned box (2D or 3D; 2D boxes keep z = 0).
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Aabb {
    pub lo: V3,
    pub hi: V3,
}

impl Aabb {
    pub const EMPTY: Aabb = Aabb { lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3] };

    pub fn new(lo: V3, hi: V3) -> Aabb {
        Aabb { lo, hi }
    }
    pub fn is_empty(&self) -> bool {
        (0..3).any(|k| self.lo[k] > self.hi[k])
    }
    pub fn grow(&mut self, p: V3) {
        for (k, &x) in p.iter().enumerate() {
            self.lo[k] = self.lo[k].min(x);
            self.hi[k] = self.hi[k].max(x);
        }
    }
    pub fn union(&self, o: &Aabb) -> Aabb {
        let mut r = *self;
        if !o.is_empty() {
            r.grow(o.lo);
            r.grow(o.hi);
        }
        r
    }
    pub fn pad(&self, d: f64) -> Aabb {
        Aabb { lo: sub3(self.lo, [d; 3]), hi: add3(self.hi, [d; 3]) }
    }
    pub fn overlaps(&self, o: &Aabb) -> bool {
        (0..3).all(|k| self.lo[k] <= o.hi[k] && o.lo[k] <= self.hi[k])
    }
    pub fn contains(&self, p: V3) -> bool {
        (0..3).all(|k| p[k] >= self.lo[k] && p[k] <= self.hi[k])
    }
    pub fn diagonal(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            norm3(sub3(self.hi, self.lo))
        }
    }
    pub fn center(&self) -> V3 {
        mul3(add3(self.lo, self.hi), 0.5)
    }
}

/// Similarity transform of the plane: `p -> m p + t` with `m = s * R` or `s * R * mirror`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Xf2 {
    pub m: [[f64; 2]; 2],
    pub t: V2,
}

impl Xf2 {
    pub const ID: Xf2 = Xf2 { m: [[1.0, 0.0], [0.0, 1.0]], t: [0.0; 2] };

    pub fn apply(&self, p: V2) -> V2 {
        [self.m[0][0] * p[0] + self.m[0][1] * p[1] + self.t[0], self.m[1][0] * p[0] + self.m[1][1] * p[1] + self.t[1]]
    }
    /// Linear part only (directions).
    pub fn apply_dir(&self, p: V2) -> V2 {
        [self.m[0][0] * p[0] + self.m[0][1] * p[1], self.m[1][0] * p[0] + self.m[1][1] * p[1]]
    }
    pub fn det(&self) -> f64 {
        self.m[0][0] * self.m[1][1] - self.m[0][1] * self.m[1][0]
    }
    /// Uniform scale factor.
    pub fn scale(&self) -> f64 {
        self.det().abs().sqrt()
    }
    pub fn mirrored(&self) -> bool {
        self.det() < 0.0
    }
    /// `self` applied after `o`.
    pub fn after(&self, o: &Xf2) -> Xf2 {
        let (a, b) = (self.m, o.m);
        let m = [
            [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
            [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
        ];
        Xf2 { m, t: self.apply(o.t) }
    }
    pub fn inverse(&self) -> Xf2 {
        let d = self.det();
        let m = [[self.m[1][1] / d, -self.m[0][1] / d], [-self.m[1][0] / d, self.m[0][0] / d]];
        let inv = Xf2 { m, t: [0.0; 2] };
        Xf2 { m, t: mul2(inv.apply(self.t), -1.0) }
    }
    pub fn translate(t: V2) -> Xf2 {
        Xf2 { m: Xf2::ID.m, t }
    }
    pub fn rotate(rad: f64) -> Xf2 {
        let (s, c) = rad.sin_cos();
        Xf2 { m: [[c, -s], [s, c]], t: [0.0; 2] }
    }
    pub fn scaling(s: f64) -> Xf2 {
        Xf2 { m: [[s, 0.0], [0.0, s]], t: [0.0; 2] }
    }
    /// Reflection across the line `p[axis] = at`.
    pub fn mirror(axis: usize, at: f64) -> Xf2 {
        let mut m = Xf2::ID.m;
        m[axis][axis] = -1.0;
        let mut t = [0.0; 2];
        t[axis] = 2.0 * at;
        Xf2 { m, t }
    }
}

/// 64-bit content hasher (splitmix-style mixing); deterministic across runs and platforms.
#[derive(Clone, Copy, Debug)]
pub struct Hasher64(pub u64);

impl Hasher64 {
    pub fn new() -> Self {
        Hasher64(0x9E37_79B9_7F4A_7C15)
    }
    pub fn u64(&mut self, v: u64) -> &mut Self {
        let mut z = (self.0 ^ v).wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        self.0 = z ^ (z >> 31);
        self
    }
    pub fn f64(&mut self, v: f64) -> &mut Self {
        self.u64(if v == 0.0 { 0 } else { v.to_bits() })
    }
    pub fn str(&mut self, s: &str) -> &mut Self {
        self.u64(s.len() as u64);
        for b in s.bytes() {
            self.u64(b as u64);
        }
        self
    }
    pub fn finish(&self) -> u64 {
        self.0
    }
}

impl Default for Hasher64 {
    fn default() -> Self {
        Self::new()
    }
}

/// Placement of a solid: a plane similarity plus an affine map of z (`z -> zs * z + zt`, `zs = +-scale`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Xf3 {
    pub xy: Xf2,
    pub zs: f64,
    pub zt: f64,
}

impl Xf3 {
    pub const ID: Xf3 = Xf3 { xy: Xf2::ID, zs: 1.0, zt: 0.0 };
    /// `self` applied after `o`.
    pub fn after(&self, o: &Xf3) -> Xf3 {
        Xf3 { xy: self.xy.after(&o.xy), zs: self.zs * o.zs, zt: self.zs * o.zt + self.zt }
    }
    pub fn interval(&self, z: (f64, f64)) -> (f64, f64) {
        let (a, b) = (self.zs * z.0 + self.zt, self.zs * z.1 + self.zt);
        (a.min(b), a.max(b))
    }
    pub fn key(&self) -> [u64; 8] {
        let f = |v: f64| if v == 0.0 { 0 } else { v.to_bits() };
        [
            self.xy.m[0][0],
            self.xy.m[0][1],
            self.xy.m[1][0],
            self.xy.m[1][1],
            self.xy.t[0],
            self.xy.t[1],
            self.zs,
            self.zt,
        ]
        .map(f)
    }
}

/// General affine map of space (rotations about any axis, mirrors, uniform scale, translation), used by polyhedral solids.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aff3 {
    pub m: [[f64; 3]; 3],
    pub t: V3,
}

impl Aff3 {
    pub const ID: Aff3 = Aff3 { m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], t: [0.0; 3] };
    pub fn apply(&self, p: V3) -> V3 {
        let r = |i: usize| self.m[i][0] * p[0] + self.m[i][1] * p[1] + self.m[i][2] * p[2] + self.t[i];
        [r(0), r(1), r(2)]
    }
    pub fn det(&self) -> f64 {
        let m = &self.m;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
    /// `self` applied after `o`.
    pub fn after(&self, o: &Aff3) -> Aff3 {
        let mut m = [[0.0; 3]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, x) in row.iter_mut().enumerate() {
                *x = (0..3).map(|k| self.m[i][k] * o.m[k][j]).sum();
            }
        }
        Aff3 { m, t: self.apply(o.t) }
    }
    pub fn translate(t: V3) -> Aff3 {
        Aff3 { t, ..Aff3::ID }
    }
    pub fn scaling(s: f64, about: V3) -> Aff3 {
        Aff3 { m: [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]], t: mul3(about, 1.0 - s) }
    }
    /// Reflection across the plane `p[axis] = at`.
    pub fn mirror(axis: usize, at: f64) -> Aff3 {
        let mut a = Aff3::ID;
        a.m[axis][axis] = -1.0;
        a.t[axis] = 2.0 * at;
        a
    }
    /// Rotation by `rad` about the line through `about` along the unit vector `k` (Rodrigues).
    pub fn rotation(k: V3, rad: f64, about: V3) -> Aff3 {
        let (s, c) = rad.sin_cos();
        let m = [
            [c + k[0] * k[0] * (1.0 - c), k[0] * k[1] * (1.0 - c) - k[2] * s, k[0] * k[2] * (1.0 - c) + k[1] * s],
            [k[1] * k[0] * (1.0 - c) + k[2] * s, c + k[1] * k[1] * (1.0 - c), k[1] * k[2] * (1.0 - c) - k[0] * s],
            [k[2] * k[0] * (1.0 - c) - k[1] * s, k[2] * k[1] * (1.0 - c) + k[0] * s, c + k[2] * k[2] * (1.0 - c)],
        ];
        let r = Aff3 { m, t: [0.0; 3] };
        Aff3::translate(about).after(&r).after(&Aff3::translate(mul3(about, -1.0)))
    }
    pub fn key(&self) -> [u64; 12] {
        let f = |v: f64| if v == 0.0 { 0 } else { v.to_bits() };
        [
            self.m[0][0],
            self.m[0][1],
            self.m[0][2],
            self.m[1][0],
            self.m[1][1],
            self.m[1][2],
            self.m[2][0],
            self.m[2][1],
            self.m[2][2],
            self.t[0],
            self.t[1],
            self.t[2],
        ]
        .map(f)
    }
}
