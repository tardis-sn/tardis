

set -ex



pytest -k "not (_not_a_real_test or (test_imageshow and test_show) or test_strip_ycbcr_jpeg_2x2_sampling or test_tiff_crashes[Tests/images/crash-81154a65438ba5aaeca73fd502fa4850fbde60f8.tif])" --durations=50
exit 0
