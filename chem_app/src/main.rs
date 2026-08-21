#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_title("Chem Workbench"),
        //renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };

    eframe::run_native(
        "Chem Workbench",
        options,
        Box::new(|cc| {
            //cc.egui_ctx.set_pixels_per_point(1.25);
            Ok(Box::new(chem_app::WorkbenchApp::new(cc)))
        }),
    )
}

#[cfg(target_arch = "wasm32")]
fn main() {
    use wasm_bindgen::JsCast as _;

    eframe::WebLogger::init(log::LevelFilter::Info).ok();
    console_error_panic_hook::set_once();

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async {
        let document = web_sys::window()
            .expect("no window")
            .document()
            .expect("no document");

        let canvas = document
            .get_element_by_id(chem_app::CANVAS_ID)
            .expect("failed to find canvas with id `the_canvas_id`")
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .expect("element with id `the_canvas_id` is not a canvas");

        let start_result = eframe::WebRunner::new()
            .start(
                canvas,
                web_options,
                Box::new(|cc| Ok(Box::new(chem_app::WorkbenchApp::new(cc)))),
            )
            .await;

        let loading_overlay = document.get_element_by_id("loading_overlay");
        match start_result {
            // The app is about to render its first real frame — the plain
            // HTML/CSS loading overlay (shown since page load, well before
            // the wasm bundle finished downloading) has done its job.
            Ok(()) => {
                if let Some(el) = loading_overlay {
                    el.remove();
                }
            }
            // Leave the overlay up but show why, rather than a silently
            // blank page.
            Err(e) => {
                log::error!("Failed to start eframe: {e:?}");
                if let Some(el) = loading_overlay {
                    el.set_text_content(Some(&format!("Failed to start Chem Workbench: {e:?}")));
                }
            }
        }
    });
}
