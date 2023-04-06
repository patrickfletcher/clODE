licenses(["notice"])  # Apache 2

#cc_library(
#    name = "python",
#    hdrs = glob([
#        "Include/*.h",
#        "PC/*.h",
#        "Modules/**/*.h",
#    ]),
#    includes = [
#        "Include",
#        "PC",
#    ],
#    visibility = ["//visibility:public"],
#)

#cc_library(
#    name = "python",
#    #    hdrs = glob([
#    #        "Include/*.h",
#    #        "PC/*.h",
#    #        "Modules/**/*.h",
#    #    ]),
#    includes = glob(
#        ["Include", "PC", "Modules", "Modules/*/*"],
#        exclude_directories = 0,
#        exclude = ["**/*.c", "**/*.h"],
#    ),
#    deps = ["python-internal"],
#    visibility = ["//visibility:public"],
#)
#
#cc_import(
#    name = "python-internal",
#    hdrs = glob([
#        "Include/*.h",
#        "PC/*.h",
#        "Modules/**/*.h",
#    ]),
#    #     static_library = ...,
#    visibility = ["//visibility:private"],
#)

cc_library(
    name = "python",
    #    hdrs = glob([
    #        "Include/*.h",
    #        "PC/*.h",
    #        "Modules/**/*.h",
    #    ]),
    includes =
        #glob(["Include", "PC", "Modules"]),
        #glob(["Include", "PC", "Modules", "Modules/**/*"]),
        ["Include", "PC", "Modules", "Modules/_decimal/libmpdec", "Include/cpython"],
    #["/Library/Framework/Python.framework/Versions/Current/Headers/"],
    deps = ["python-internal"],
    visibility = ["//visibility:public"],
)

cc_import(
    name = "python-internal",
    hdrs = glob([
        "Include/*.h",
        "PC/*.h",
        "Modules/**/*.h",
    ]),
    #     static_library = ...,
    visibility = ["//visibility:private"],
)
