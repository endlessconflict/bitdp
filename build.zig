const std = @import("std");

pub fn build(b: *std.Build) void {
    const query = b.standardTargetOptionsQueryOnly(.{});
    const target = b.resolveTargetQuery(query);
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("bitdp", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const tests = b.addTest(.{ .root_module = mod });
    b.step("test", "Run tests").dependOn(&b.addRunArtifact(tests).step);

    const mve = b.addExecutable(.{
        .name = "bitdp-mve",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/mve.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "bitdp", .module = mod }},
        }),
    });
    const run = b.addRunArtifact(mve);
    if (b.args) |args| run.addArgs(args);
    b.step("mve", "Verify and benchmark derived kernels (use -Doptimize=ReleaseFast)").dependOn(&run.step);

    const bench = b.addExecutable(.{
        .name = "bitdp-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "bitdp", .module = mod }},
        }),
    });
    b.step("bench", "Build the benchmark driver (see bench/README.md)").dependOn(&b.addInstallArtifact(bench, .{}).step);

    const scan = b.addExecutable(.{
        .name = "bitdp-scan",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/scan.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "bitdp", .module = mod }},
        }),
    });
    b.step("scan", "Build the genome scan demo (see bench/README.md)").dependOn(&b.addInstallArtifact(scan, .{}).step);

    const opts = b.addOptions();
    opts.addOption(u32, "scheme", b.option(u32, "scheme", "Ablation: scheme index in bench/ablation.zig") orelse 0);
    const ablation = b.addExecutable(.{
        .name = "bitdp-ablation",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/ablation.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "bitdp", .module = mod },
                .{ .name = "options", .module = opts.createModule() },
            },
        }),
    });
    b.step("ablation", "Print op counts per builder configuration for -Dscheme=N").dependOn(&b.addRunArtifact(ablation).step);

    // C interface: libbitdp (static and shared) with include/bitdp.h.
    // The library is meant to be copied to other machines, so on x86-64 it targets
    // x86_64_v3 (AVX2) unless -Dcpu is given; everything else builds for the native CPU.
    var lib_query = query;
    if (query.cpu_model == .determined_by_arch_os and (query.cpu_arch orelse @import("builtin").cpu.arch) == .x86_64)
        lib_query.cpu_model = .{ .explicit = &std.Target.x86.cpu.x86_64_v3 };
    const lib_target = b.resolveTargetQuery(lib_query);
    // PIC so the static library links into the position-independent executables C toolchains build by default.
    const capi = b.createModule(.{ .root_source_file = b.path("src/c_api.zig"), .target = lib_target, .optimize = optimize, .pic = true });
    const lib_static = b.addLibrary(.{ .name = "bitdp", .root_module = capi, .linkage = .static });
    const lib_shared = b.addLibrary(.{
        .name = "bitdp",
        .root_module = b.createModule(.{ .root_source_file = b.path("src/c_api.zig"), .target = lib_target, .optimize = optimize }),
        .linkage = .dynamic,
    });
    const lib = b.step("lib", "Build libbitdp (static and shared) and install include/bitdp.h");
    lib.dependOn(&b.addInstallArtifact(lib_static, .{}).step);
    lib.dependOn(&b.addInstallArtifact(lib_shared, .{ .implib_dir = .disabled }).step);
    lib.dependOn(&b.addInstallFileWithDir(b.path("include/bitdp.h"), .header, "bitdp.h").step);

    // The C test, built once as C99 and once as C++, linked against the static library.
    const test_c = b.step("test-c", "Run the C interface tests (as C and as C++)");
    for ([_]std.Build.Module.CSourceLanguage{ .c, .cpp }) |lang| {
        const m = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true, .link_libcpp = lang == .cpp });
        m.addCSourceFile(.{ .file = b.path("tests/c/test_bitdp.c"), .flags = if (lang == .c) &.{ "-std=c99", "-Wall", "-Wextra", "-Werror" } else &.{ "-Wall", "-Wextra", "-Werror" }, .language = lang });
        m.addIncludePath(b.path("include"));
        m.linkLibrary(lib_static);
        const exe = b.addExecutable(.{ .name = if (lang == .c) "test-c" else "test-cpp", .root_module = m });
        test_c.dependOn(&b.addRunArtifact(exe).step);
    }
}
