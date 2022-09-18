load(
    "@bazel_tools//tools/cpp:lib_cc_configure.bzl",
    "escape_string",
    "which",
)
load("@bazel_tools//tools/build_defs/cc:action_names.bzl", "ACTION_NAMES")
load(
    "@bazel_tools//tools/cpp:cc_toolchain_config_lib.bzl",
    "action_config",
    "feature",
    "flag_group",
    "flag_set",
    "tool",
    "tool_path",
    "with_feature_set",
)

def _find_tool_in_applications(repository_ctx, tool, overriden_tools):
    """Find a tool for repository, taking overriden tools into account."""
    if tool in overriden_tools:
        return overriden_tools[tool]
    return which(repository_ctx, tool, "/Applications/" + tool)

def _get_tool_paths(repository_ctx, overriden_tools):
    """Compute the %-escaped path to the various tools"""
    return dict({
        k: escape_string(_find_tool_in_applications(repository_ctx, k, overriden_tools))
        for k in [
            "mex",
        ]
    }.items())

all_link_actions = [
    # NEW
    ACTION_NAMES.cpp_link_executable,
    ACTION_NAMES.cpp_link_dynamic_library,
    ACTION_NAMES.cpp_link_static_library,
    ACTION_NAMES.cpp_link_nodeps_dynamic_library,
]

#
all_compile_actions = [
    ACTION_NAMES.assemble,
    ACTION_NAMES.preprocess_assemble,
    ACTION_NAMES.c_compile,
    ACTION_NAMES.cpp_compile,
    ACTION_NAMES.cpp_module_compile,
    ACTION_NAMES.objc_compile,
    ACTION_NAMES.objcpp_compile,
    ACTION_NAMES.cpp_header_parsing,
    ACTION_NAMES.clif_match,
]

def _impl(ctx):
    print("ctx.attr", ctx.attr)
    print("ctx.bin_dir", ctx.bin_dir.path)
    print("ctx.build_file_path", ctx.build_file_path)
    print("ctx.genfiles_dir", ctx.genfiles_dir.path)
    print(dir(ctx.label))
    output_location = "//" + ctx.label.package + ":" + ctx.label.name
    print("Ouput location", output_location)
    print("ctx.expand_location", ctx.expand_location("Test %(" + output_location + ")"))
    tool_paths = [
        tool_path(
            name = "gcc",
            #path = "/Applications/MATLAB_R2022a.app/bin/mex",
            path = "mex_osx.sh",
        ),
        tool_path(
            name = "ld",
            path = "/bin/true",
        ),
        tool_path(
            name = "ar",
            path = "/bin/false",
        ),
        tool_path(
            name = "cpp",
            path = "/bin/false",
        ),
        tool_path(
            name = "gcov",
            path = "/bin/false",
        ),
        tool_path(
            name = "nm",
            path = "/bin/false",
        ),
        tool_path(
            name = "objdump",
            path = "/bin/false",
        ),
        tool_path(
            name = "strip",
            path = "/bin/false",
        ),
    ]

    compile_action_configs = [
        action_config(
            action_name = action_name,
            tools = [
                tool(
                    path = "mex_osx.sh",
                ),
            ],
        )
        for action_name in all_compile_actions
    ]
    link_action_configs = [
        action_config(
            action_name = action_name,
            tools = [
                tool(
                    #path = "/Applications/MATLAB_R2022a.app/bin/mex",
                    path = "/usr/bin/ld",
                ),
            ],
        )
        for action_name in all_link_actions
    ]
    action_configs = [
        action_config(
            action_name = ACTION_NAMES.cpp_compile,
            tools = [
                tool(
                    path = "mex_osx.sh",
                ),
            ],
        ),
        action_config(
            action_name = ACTION_NAMES.c_compile,
            tools = [
                tool(
                    path = "mex_osx.sh",
                ),
            ],
        ),
        action_config(
            action_name = ACTION_NAMES.cpp_module_compile,
            tools = [
                tool(
                    path = "mex_osx.sh",
                ),
            ],
        ),
    ]

    features = [
        feature(
            name = "no_legacy_features",
            enabled = True,
        ),
        feature(
            name = "includes",
            enabled = True,
            flag_sets = [
                flag_set(
                    actions = all_compile_actions,
                    flag_groups = [
                        flag_group(
                            iterate_over = "include_paths",
                            flags = ["-I%{include_paths}"],
                        ),
                        flag_group(
                            iterate_over = "quote_include_paths",
                            flags = ["-I%{quote_include_paths}"],
                        ),
                        flag_group(
                            iterate_over = "system_include_paths",
                            flags = ["-I%{system_include_paths}"],
                        ),
                    ],
                ),
            ],
        ),
        feature(
            name = "source_files",
            enabled = True,
            flag_sets = [
                flag_set(
                    actions = all_compile_actions,
                    flag_groups = [
                        flag_group(
                            flags = ["%{source_file}"],
                        ),
                    ],
                ),
            ],
        ),
        feature(
            name = "output_object_files",
            enabled = True,
            flag_sets = [
                flag_set(
                    actions = all_compile_actions,
                    flag_groups = [
                        flag_group(
                            flags = ["-c", "-v"],
                        ),
                        flag_group(
                            expand_if_available = "output_file",
                            flags = ["-o", "%{output_file}"],
                        ),
                    ],
                ),
            ],
        ),
        feature(
            name = "cxx_flags",
            enabled = True,
            flag_sets = [
                flag_set(
                    actions = all_compile_actions,
                    flag_groups = ([
                        flag_group(
                            flags = [
                                "CXXFLAGS=$CXXFLAGS -MD -MF %{dependency_file}",
                            ],
                            expand_if_available = "dependency_file",
                        ),
                    ]),
                ),
            ],
        ),
        feature(
            name = "dependency_file",
            enabled = True,
            flag_sets = [
                flag_set(
                    actions = all_compile_actions,
                    flag_groups = ([
                        flag_group(
                            flags = [
                                "-dependency_file",
                                "%{dependency_file}",
                            ],
                            expand_if_available = "dependency_file",
                        ),
                    ]),
                ),
            ],
        ),
        feature(
            name = "linker_param_file",
            enabled = True,
            flag_sets = [
                flag_set(
                    actions = all_link_actions + [
                        ACTION_NAMES.cpp_link_static_library,
                        ACTION_NAMES.objc_archive,
                        ACTION_NAMES.objc_fully_link,
                        ACTION_NAMES.objc_executable,
                        ACTION_NAMES.objcpp_executable,
                    ],
                    flag_groups = [
                        flag_group(
                            #flags = ["@%{linker_param_file}"],
                            flags = ["bazel-out/mez-fastbuild/bin/example/libtest_class.a-2.params"],
                            #expand_if_available = "linker_param_file",
                        ),
                    ],
                ),
            ],
        ),
    ]

    return cc_common.create_cc_toolchain_config_info(
        ctx = ctx,
        features = features,
        action_configs = compile_action_configs + link_action_configs,
        toolchain_identifier = "mex-toolchain",
        host_system_name = "local",
        target_system_name = "local",
        target_cpu = "mex",
        target_libc = "unknown",
        compiler = "mex_osx.sh",
        abi_version = "unknown",
        abi_libc_version = "unknown",
        tool_paths = tool_paths,
    )

cc_toolchain_config = rule(
    implementation = _impl,
    attrs = {},
    provides = [CcToolchainConfigInfo],
)
