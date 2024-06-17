using Aqua

Aqua.test_all(
    StereoSSAM; ambiguities=(recursive = false), deps_compat=(check_weakdeps = false)
)
