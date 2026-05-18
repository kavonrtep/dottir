use dottir_core::{pick_auto_zoom, snap_zoom_to_period_divisor};

#[test]
fn plain_auto_zoom_is_unchanged() {
    assert_eq!(pick_auto_zoom(27_000, 27_000, 4096), 7);
}

#[test]
fn snap_keeps_exact_divisor() {
    assert_eq!(snap_zoom_to_period_divisor(11, &[451], 2.0), 11);
}

#[test]
fn snap_prefers_nearby_coarser_divisor() {
    // 451 = 11 * 41. For a base auto-fit zoom of 22, choose 41 rather
    // than 11 so auto-fit does not silently allocate a larger pixelmap.
    assert_eq!(snap_zoom_to_period_divisor(22, &[451], 2.0), 41);
}

#[test]
fn snap_uses_common_divisor_for_multiple_periods() {
    assert_eq!(snap_zoom_to_period_divisor(20, &[450, 900], 2.0), 25);
}

#[test]
fn snap_falls_back_without_reasonable_period() {
    assert_eq!(snap_zoom_to_period_divisor(22, &[], 2.0), 22);
    assert_eq!(snap_zoom_to_period_divisor(22, &[451], 1.1), 22);
    assert_eq!(snap_zoom_to_period_divisor(22, &[1], 2.0), 22);
}
