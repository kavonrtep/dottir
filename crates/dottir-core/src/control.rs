//! Cancellation and progress hooks shared by the exact and large engines.

use std::sync::atomic::{AtomicBool, Ordering};

use crate::DottirError;

/// Cancellation and progress hooks. Progress is reported as a fraction in
/// `[0, 1]` from worker threads.
#[derive(Default, Clone, Copy)]
pub struct Control<'a> {
    pub cancel: Option<&'a AtomicBool>,
    pub progress: Option<&'a (dyn Fn(f32) + Sync)>,
}

impl std::fmt::Debug for Control<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Control")
            .field("cancel", &self.cancel.map(|c| c.load(Ordering::Relaxed)))
            .field("progress", &self.progress.is_some())
            .finish()
    }
}

impl Control<'_> {
    pub(crate) fn cancelled(&self) -> bool {
        self.cancel.is_some_and(|c| c.load(Ordering::Relaxed))
    }
    pub(crate) fn check(&self) -> Result<(), DottirError> {
        if self.cancelled() {
            Err(DottirError::Cancelled)
        } else {
            Ok(())
        }
    }
    pub(crate) fn report(&self, f: f32) {
        if let Some(p) = self.progress {
            p(f.clamp(0.0, 1.0));
        }
    }
}
