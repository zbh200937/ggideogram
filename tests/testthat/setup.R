# Prevent bare ggplot2::ggplotGrob() calls in testthat from writing Rplots.pdf
# to disk; tests that explicitly assert no device is opened override this option
# locally with a throwing device callback.
old_ideogram_test_device <- options(
  device = function(...) grDevices::pdf(file = NULL)
)
teardown_target <- tryCatch(testthat::teardown_env(), error = function(e) NULL)
if (!is.null(teardown_target)) {
  reg.finalizer(
    teardown_target,
    function(...) {
      while (grDevices::dev.cur() > 1L) grDevices::dev.off()
      options(old_ideogram_test_device)
    },
    onexit = TRUE
  )
}
