//! Deferred work whose result is picked up from the UI loop.
//!
//! Three operations — dataset fingerprinting, query fingerprinting and
//! similarity search — can run on the GPU, which means they're async on
//! wasm32: `eframe`'s `update()` is synchronous and can't await, so the work is
//! spawned and its result collected on a later frame. Each had its own
//! `cfg`-split dispatch pair, its own `Rc<RefCell<Option<..>>>` slot, and its
//! own polling block in `update()`, all following the same shape.
//!
//! [`Task`] holds that shape once. What it deliberately doesn't try to hold is
//! the operations themselves: their inputs and outputs have nothing in common
//! (a slice of molecules to a vector of fingerprints, one molecule to one
//! fingerprint, a query plus targets to ranked hits), and what each result
//! *means* to the app differs too. The repetition was in the plumbing, not the
//! work.

use std::future::Future;
use web_time::Instant;

#[cfg(target_arch = "wasm32")]
use std::{cell::RefCell, rc::Rc};

/// A result and how long it took, in milliseconds.
pub type Timed<T> = (anyhow::Result<T>, f64);

/// Work that may finish later, with its result collected by [`Task::poll`].
///
/// On wasm32 the work is spawned and lands in a shared slot; on native it runs
/// to completion immediately and waits in the same place. Either way the caller
/// starts it in one place and handles the result in another, so a single code
/// path serves both platforms.
pub struct Task<T> {
    #[cfg(target_arch = "wasm32")]
    slot: Rc<RefCell<Option<Timed<T>>>>,
    #[cfg(not(target_arch = "wasm32"))]
    ready: Option<Timed<T>>,
}

impl<T: 'static> Default for Task<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: 'static> Task<T> {
    pub fn new() -> Self {
        Self {
            #[cfg(target_arch = "wasm32")]
            slot: Rc::new(RefCell::new(None)),
            #[cfg(not(target_arch = "wasm32"))]
            ready: None,
        }
    }

    /// Starts `future`, timing it.
    ///
    /// The future has to own everything it touches (`'static`): on wasm32 it
    /// outlives this call, so it can't borrow from the app. In practice callers
    /// clone what they need first — cheaply, since the search engine's GPU
    /// handles are `Arc`-backed.
    ///
    /// Starting a second task before the first is polled discards the first
    /// result. That's what should happen: it belongs to a superseded request,
    /// and applying it would overwrite the newer one.
    pub fn start<F>(&mut self, future: F)
    where
        F: Future<Output = anyhow::Result<T>> + 'static,
    {
        #[cfg(target_arch = "wasm32")]
        {
            let slot = self.slot.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let started = Instant::now();
                let result = future.await;
                let elapsed_ms = started.elapsed().as_secs_f64() * 1000.0;
                *slot.borrow_mut() = Some((result, elapsed_ms));
            });
        }

        #[cfg(not(target_arch = "wasm32"))]
        {
            // Native has threads, so blocking here only stalls this one — and
            // it's what the UI did before. The result waits to be polled rather
            // than being applied inline, so both platforms follow the same
            // path and a bug in the collection step can't hide on one of them.
            let started = Instant::now();
            let result = pollster::block_on(future);
            let elapsed_ms = started.elapsed().as_secs_f64() * 1000.0;
            self.ready = Some((result, elapsed_ms));
        }
    }

    /// Takes the result if one is waiting, leaving the task empty.
    pub fn poll(&mut self) -> Option<Timed<T>> {
        #[cfg(target_arch = "wasm32")]
        {
            self.slot.borrow_mut().take()
        }

        #[cfg(not(target_arch = "wasm32"))]
        {
            self.ready.take()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nothing_to_poll_before_starting() {
        let mut task: Task<u32> = Task::new();
        assert!(task.poll().is_none());
    }

    #[test]
    fn test_result_is_available_after_starting() {
        let mut task = Task::new();
        task.start(async { Ok(7u32) });

        let (result, elapsed_ms) = task.poll().expect("a result should be waiting");
        assert_eq!(result.unwrap(), 7);
        assert!(elapsed_ms >= 0.0);
    }

    #[test]
    fn test_polling_consumes_the_result() {
        let mut task = Task::new();
        task.start(async { Ok(1u32) });

        assert!(task.poll().is_some());
        // Applying the same result twice would double-count it.
        assert!(task.poll().is_none());
    }

    #[test]
    fn test_errors_are_carried_through() {
        let mut task: Task<u32> = Task::new();
        task.start(async { Err(anyhow::anyhow!("boom")) });

        let (result, _) = task.poll().expect("a result should be waiting");
        assert_eq!(result.unwrap_err().to_string(), "boom");
    }

    #[test]
    fn test_restarting_discards_the_earlier_result() {
        let mut task = Task::new();
        task.start(async { Ok(1u32) });
        // The first result belongs to a superseded request; applying it would
        // overwrite the newer one.
        task.start(async { Ok(2u32) });

        let (result, _) = task.poll().expect("a result should be waiting");
        assert_eq!(result.unwrap(), 2);
        assert!(task.poll().is_none());
    }
}
