//! Small lossless byte-stream compressor: an adaptive binary range coder (the carry-propagating
//! scheme used by LZMA) driven by an order-0 or order-1 byte model. Pure Rust, no dependencies.
//!
//! It is the last stage for every array in a `.fmsh` file. Each stream is stored with whichever of
//! raw / order-0 / order-1 is smallest, so incompressible data costs only a few bytes extra.

const PROB_BITS: u32 = 15;
const PROB_INIT: u16 = 1 << (PROB_BITS - 1);
const MOVE_BITS: u32 = 5;
const TOP: u32 = 1 << 24;

struct Encoder {
    low: u64,
    range: u32,
    cache: u8,
    cache_size: u64,
    out: Vec<u8>,
}

impl Encoder {
    fn new() -> Self {
        Encoder { low: 0, range: 0xFFFF_FFFF, cache: 0, cache_size: 1, out: Vec::new() }
    }

    fn shift_low(&mut self) {
        if (self.low as u32) < 0xFF00_0000 || (self.low >> 32) != 0 {
            let carry = (self.low >> 32) as u8;
            let mut temp = self.cache;
            loop {
                self.out.push(temp.wrapping_add(carry));
                temp = 0xFF;
                self.cache_size -= 1;
                if self.cache_size == 0 {
                    break;
                }
            }
            self.cache = ((self.low >> 24) & 0xFF) as u8;
        }
        self.cache_size += 1;
        self.low = (self.low & 0x00FF_FFFF) << 8;
    }

    fn bit(&mut self, p: &mut u16, bit: u32) {
        let bound = (self.range >> PROB_BITS) * *p as u32;
        if bit == 0 {
            self.range = bound;
            *p += ((1 << PROB_BITS) - *p) >> MOVE_BITS;
        } else {
            self.low += bound as u64;
            self.range -= bound;
            *p -= *p >> MOVE_BITS;
        }
        while self.range < TOP {
            self.range <<= 8;
            self.shift_low();
        }
    }

    fn finish(mut self) -> Vec<u8> {
        for _ in 0..5 {
            self.shift_low();
        }
        self.out
    }
}

struct Decoder<'a> {
    b: &'a [u8],
    i: usize,
    range: u32,
    code: u32,
}

impl<'a> Decoder<'a> {
    fn new(b: &'a [u8]) -> Self {
        let mut d = Decoder { b, i: 1, range: 0xFFFF_FFFF, code: 0 };
        for _ in 0..4 {
            d.code = (d.code << 8) | d.next() as u32;
        }
        d
    }

    fn next(&mut self) -> u8 {
        let v = self.b.get(self.i).copied().unwrap_or(0);
        self.i += 1;
        v
    }

    fn bit(&mut self, p: &mut u16) -> u32 {
        let bound = (self.range >> PROB_BITS) * *p as u32;
        let bit = if self.code < bound {
            self.range = bound;
            *p += ((1 << PROB_BITS) - *p) >> MOVE_BITS;
            0
        } else {
            self.code -= bound;
            self.range -= bound;
            *p -= *p >> MOVE_BITS;
            1
        };
        while self.range < TOP {
            self.range <<= 8;
            self.code = (self.code << 8) | self.next() as u32;
        }
        bit
    }
}

/// Context model of the byte coder.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    /// One adaptive distribution for all bytes.
    Order0,
    /// One distribution per previous byte.
    Order1,
}

impl Model {
    fn contexts(self) -> usize {
        match self {
            Model::Order0 => 1,
            Model::Order1 => 256,
        }
    }
}

/// Compress `data` with the given model (the length is not stored).
pub fn compress(data: &[u8], model: Model) -> Vec<u8> {
    let mut probs = vec![PROB_INIT; model.contexts() * 256];
    let mut enc = Encoder::new();
    let mut prev = 0usize;
    for &byte in data {
        let base = if model == Model::Order1 { prev * 256 } else { 0 };
        let mut node = 1usize;
        for i in (0..8).rev() {
            let bit = ((byte >> i) & 1) as u32;
            enc.bit(&mut probs[base + node], bit);
            node = (node << 1) | bit as usize;
        }
        prev = byte as usize;
    }
    enc.finish()
}

/// Inverse of [`compress`]; `n` is the number of bytes to produce.
pub fn decompress(data: &[u8], n: usize, model: Model) -> Vec<u8> {
    let mut probs = vec![PROB_INIT; model.contexts() * 256];
    let mut dec = Decoder::new(data);
    let mut out = Vec::with_capacity(n);
    let mut prev = 0usize;
    for _ in 0..n {
        let base = if model == Model::Order1 { prev * 256 } else { 0 };
        let mut node = 1usize;
        for _ in 0..8 {
            let bit = dec.bit(&mut probs[base + node]);
            node = (node << 1) | bit as usize;
        }
        let byte = (node & 0xFF) as u8;
        out.push(byte);
        prev = byte as usize;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *seed >> 33
    }

    #[test]
    fn roundtrip_all_models() {
        let mut s = 7u64;
        let random: Vec<u8> = (0..5000).map(|_| lcg(&mut s) as u8).collect();
        let text: Vec<u8> =
            b"the quick brown fox jumps over the lazy dog ".iter().cycle().take(4000).copied().collect();
        let zeros = vec![0u8; 10_000];
        let ones = vec![0xFFu8; 10_000]; // exercises the carry path
        for data in [vec![], vec![42], random, text, zeros, ones] {
            for m in [Model::Order0, Model::Order1] {
                let c = compress(&data, m);
                assert_eq!(decompress(&c, data.len(), m), data, "{m:?} len {}", data.len());
            }
        }
    }

    #[test]
    fn skewed_data_shrinks_and_random_does_not_blow_up() {
        let zeros = vec![0u8; 100_000];
        assert!(compress(&zeros, Model::Order0).len() < 250, "{}", compress(&zeros, Model::Order0).len());
        let mut s = 3u64;
        let random: Vec<u8> = (0..100_000).map(|_| lcg(&mut s) as u8).collect();
        let c = compress(&random, Model::Order0);
        assert!(c.len() < 100_000 + 100_000 / 50, "random expanded to {}", c.len());
    }

    #[test]
    fn order1_beats_order0_on_correlated_bytes() {
        let data: Vec<u8> = (0..20_000u32).map(|i| (i % 7 * 31) as u8).collect();
        assert!(compress(&data, Model::Order1).len() < compress(&data, Model::Order0).len());
    }
}
