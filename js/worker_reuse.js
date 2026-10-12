// Saved snapshots remain owned by the page. If the solve can use a snapshot,
// omit its return copy only after the worker acknowledges matching live state.
// Clear the tracker on worker replacement/errors and before inverse design.
export function createWorkerReuse() {
  let cachedState = null;
  let pendingId = null;

  return {
    prepare(requestId, settings) {
      const usesSnapshot = settings.viscous && (!settings.advancedMode || settings.reuseSolution);
      const canOmit = !usesSnapshot || (cachedState && settings.reuseState === cachedState
        && settings.geometryKey === cachedState.geometryKey
        && settings.reusePanel);
      cachedState = null;
      pendingId = requestId;
      return { ...settings, reportReuseCache: true,
        reuseState: canOmit ? null : settings.reuseState };
    },
    accept(payload) {
      // A request queued after this solve (e.g. QDES) invalidates its ack even
      // if the solve's response has not yet reached the page.
      if (pendingId === null || payload.type || payload.requestId !== pendingId) return;
      pendingId = null;
      cachedState = payload.ok && payload.converged && payload.reuseCacheReady === true
        && payload.reuseState?.bl?.LVCONV === true ? payload.reuseState : null;
    },
    clear() {
      cachedState = null;
      pendingId = null;
    },
  };
}
