//! Handing an SVG to the user.
//!
//! Producing the document is [`chemdraw::svg::structure_to_svg`]; what stays
//! here is the part that cannot leave the application — a native file dialog on
//! one platform, a browser download on the other.

/// Hands the SVG to the user: a file on native, a download on web.
///
/// The same split as loading, and for the same reason — a browser has no
/// filesystem to write to, so the two platforms need different mechanisms for
/// what is conceptually one action.
#[cfg(not(target_arch = "wasm32"))]
pub fn save_svg(suggested_name: &str, contents: &str) {
    // Blocking is what the load path does too: on native it stalls this thread
    // only, and the dialog is modal to the user regardless.
    let picked = pollster::block_on(
        rfd::AsyncFileDialog::new()
            .set_file_name(suggested_name)
            .add_filter("SVG", &["svg"])
            .save_file(),
    );

    if let Some(file) = picked
        && let Err(e) = std::fs::write(file.path(), contents)
    {
        log::error!("Failed to write SVG: {e}");
    }
}

/// On web there is nothing to write to, so the file is handed to the browser as
/// a blob behind a synthetic download link.
#[cfg(target_arch = "wasm32")]
pub fn save_svg(suggested_name: &str, contents: &str) {
    use wasm_bindgen::JsCast as _;

    let Some(document) = web_sys::window().and_then(|w| w.document()) else {
        log::error!("No document to download through");
        return;
    };

    // A one-element array of one string: what Blob's constructor takes.
    let parts = js_sys::Array::new();
    parts.push(&wasm_bindgen::JsValue::from_str(contents));
    let properties = web_sys::BlobPropertyBag::new();
    properties.set_type("image/svg+xml");

    let Ok(blob) = web_sys::Blob::new_with_str_sequence_and_options(&parts, &properties) else {
        log::error!("Failed to build the SVG blob");
        return;
    };
    let Ok(url) = web_sys::Url::create_object_url_with_blob(&blob) else {
        log::error!("Failed to create a download URL");
        return;
    };

    let anchor = document
        .create_element("a")
        .ok()
        .and_then(|e| e.dyn_into::<web_sys::HtmlAnchorElement>().ok());
    if let Some(anchor) = anchor {
        anchor.set_href(&url);
        anchor.set_download(suggested_name);
        anchor.click();
    }

    // The blob stays alive as long as its URL does, so release it now the click
    // has been dispatched.
    let _ = web_sys::Url::revoke_object_url(&url);
}
