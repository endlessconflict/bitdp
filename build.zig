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
}
