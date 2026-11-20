//! Pixelmap cache bounded by memory (docs/GUI_NAVIGATION_PLAN.md N4).
//!
//! Least recently used entries are evicted first once the total size
//! exceeds the budget. An entry larger than the whole budget is not
//! cached.

use std::collections::VecDeque;

use dottir_core::large::LargePlot;
use dottir_core::DotPlot;

/// Default budget: 1 GiB.
pub const DEFAULT_BUDGET: u64 = 1 << 30;

pub struct ByteLru<K, V> {
    /// Oldest first; a hit moves the entry to the back.
    entries: VecDeque<(K, V, u64)>,
    used: u64,
    budget: u64,
}

impl<K: PartialEq, V: Clone> ByteLru<K, V> {
    pub fn new(budget: u64) -> Self {
        Self {
            entries: VecDeque::new(),
            used: 0,
            budget,
        }
    }

    pub fn get(&mut self, key: &K) -> Option<V> {
        let pos = self.entries.iter().position(|(k, _, _)| k == key)?;
        let e = self.entries.remove(pos)?;
        let v = e.1.clone();
        self.entries.push_back(e);
        Some(v)
    }

    /// Insert or replace `key`. Returns `false` (and caches nothing) when
    /// `bytes` alone exceeds the budget.
    pub fn insert(&mut self, key: K, value: V, bytes: u64) -> bool {
        self.remove(&key);
        if bytes > self.budget {
            return false;
        }
        self.used += bytes;
        self.entries.push_back((key, value, bytes));
        self.evict();
        true
    }

    pub fn set_budget(&mut self, budget: u64) {
        self.budget = budget;
        self.evict();
    }

    pub fn clear(&mut self) {
        self.entries.clear();
        self.used = 0;
    }

    pub fn used(&self) -> u64 {
        self.used
    }

    pub fn budget(&self) -> u64 {
        self.budget
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    fn remove(&mut self, key: &K) {
        if let Some(pos) = self.entries.iter().position(|(k, _, _)| k == key) {
            if let Some((_, _, b)) = self.entries.remove(pos) {
                self.used -= b;
            }
        }
    }

    fn evict(&mut self) {
        while self.used > self.budget {
            match self.entries.pop_front() {
                Some((_, _, b)) => self.used -= b,
                None => break,
            }
        }
    }
}

/// Heap bytes of a computed view: the pixelmap(s) and, for the large
/// engine, every class channel, coverage, depth bins, segments and chains.
pub fn view_bytes(plot: &DotPlot, large: Option<&LargePlot>) -> u64 {
    let px = plot.pixels.len() + plot.reverse_pixels.as_ref().map_or(0, Vec::len);
    let lp = large.map_or(0, |lp| {
        let vecs = |v: &[Vec<u8>]| v.iter().map(Vec::len).sum::<usize>();
        vecs(&lp.channels)
            + vecs(&lp.coverage)
            + (lp.q_capped.len() + lp.s_capped.len()) * 4
            + (lp.q_depth_bins.len() + lp.s_depth_bins.len())
                * std::mem::size_of::<dottir_core::large::DepthBin>()
            + lp.segments.len() * std::mem::size_of::<dottir_core::large::Segment>()
            + lp.chains.len() * std::mem::size_of::<dottir_core::large::Chain>()
    });
    (px + lp) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evicts_least_recently_used_by_bytes() {
        let mut c = ByteLru::new(100);
        assert!(c.insert('a', 1, 40));
        assert!(c.insert('b', 2, 40));
        assert_eq!(c.get(&'a'), Some(1)); // a is now the newest
        assert!(c.insert('c', 3, 40)); // 120 > 100: b goes
        assert_eq!(c.get(&'b'), None);
        assert_eq!(c.get(&'a'), Some(1));
        assert_eq!(c.get(&'c'), Some(3));
        assert_eq!(c.used(), 80);
        assert_eq!(c.len(), 2);
    }

    #[test]
    fn oversize_entry_is_not_cached() {
        let mut c = ByteLru::new(100);
        assert!(c.insert('a', 1, 60));
        assert!(!c.insert('b', 2, 101));
        assert_eq!(c.get(&'b'), None);
        assert_eq!(c.get(&'a'), Some(1), "existing entries stay");
    }

    #[test]
    fn replacing_a_key_frees_its_bytes() {
        let mut c = ByteLru::new(100);
        c.insert('a', 1, 60);
        c.insert('a', 2, 30);
        assert_eq!(c.used(), 30);
        assert_eq!(c.get(&'a'), Some(2));
    }

    #[test]
    fn shrinking_the_budget_evicts() {
        let mut c = ByteLru::new(100);
        c.insert('a', 1, 50);
        c.insert('b', 2, 50);
        c.set_budget(60);
        assert_eq!(c.len(), 1);
        assert_eq!(c.get(&'b'), Some(2));
        c.clear();
        assert_eq!((c.len(), c.used()), (0, 0));
    }
}
