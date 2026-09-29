const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
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
}
