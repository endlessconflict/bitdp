//! Minimal viable experiment: derive kernels, print them, verify 10^6 random
//! pairs against the scalar oracle, and time them against Myers' kernel.

const std = @import("std");
const bitdp = @import("bitdp");

fn printPlan(comptime s: bitdp.Scheme) void {
    const plan = comptime bitdp.derive(s);
    std.debug.print("\nscheme match={d} mismatch={d} gap={d}: differences {any}, {d} ops/word, {d} carry chain(s), in_pol={}\n", .{
        s.match, s.mismatch, s.gap, plan.vals[0..plan.k], plan.cost, plan.chains, plan.in_pol,
    });
    for (plan.nodes[0..plan.len], 0..) |n, i| {
        std.debug.print("  r{d:<3} = {s:<5} {d} {d}\n", .{ i, @tagName(n.op), n.a, n.b });
    }
    std.debug.print("  out_v {any}  out_u {any}  lvl_pol {any}\n", .{ plan.out_v[0 .. plan.k - 1], plan.out_u[0 .. plan.k - 1], plan.lvl_pol[1..plan.k] });
}

fn verify(comptime s: bitdp.Scheme, pairs: usize) !void {
    const K = bitdp.Kernel(s);
    var prng = std.Random.DefaultPrng.init(0xb17d9);
    const rnd = prng.random();
    var buf: [65]i64 = undefined;
    var p: [64]u8 = undefined;
    var t: [256]u8 = undefined;
    for (0..pairs) |_| {
        const m = rnd.intRangeAtMost(usize, 1, 64);
        const n = rnd.intRangeAtMost(usize, 0, 256);
        for (p[0..m]) |*x| x.* = "ACGT"[rnd.int(u2)];
        for (t[0..n]) |*x| x.* = "ACGT"[rnd.int(u2)];
        const want = bitdp.reference.scalar(s, p[0..m], t[0..n], &buf);
        const got = K.distance(p[0..m], t[0..n]);
        if (want != got) {
            std.debug.print("MISMATCH m={d} n={d} want={d} got={d}\n", .{ m, n, want, got });
            return error.Mismatch;
        }
    }
    std.debug.print("  verified {d} random pairs (m<=64, n<=256): bit-exact\n", .{pairs});
}

fn now(io: std.Io) i96 {
    return std.Io.Timestamp.now(io, .awake).nanoseconds;
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    printPlan(.{});
    printPlan(.{ .mismatch = 2 });
    printPlan(.{ .mismatch = 1, .gap = 2 });
    printPlan(.{ .mismatch = 3, .gap = 2 });

    try verify(.{}, 1_000_000);
    try verify(.{ .mismatch = 2 }, 200_000);
    try verify(.{ .mismatch = 1, .gap = 2 }, 200_000);
    try verify(.{ .mismatch = 3, .gap = 2 }, 200_000);

    // Throughput: m = 64 against a long random text.
    const gpa = init.gpa;
    const n: usize = 50_000_000;
    const text = try gpa.alloc(u8, n);
    defer gpa.free(text);
    var prng = std.Random.DefaultPrng.init(7);
    for (text) |*x| x.* = "ACGT"[prng.random().int(u2)];
    const pat = text[1000..1064];
    const K = bitdp.Kernel(.{});

    var t0 = now(io);
    const a = bitdp.reference.myers(pat, text);
    const t_myers = now(io) - t0;
    t0 = now(io);
    const b = K.distance(pat, text);
    const t_derived = now(io) - t0;
    std.mem.doNotOptimizeAway(a);
    std.mem.doNotOptimizeAway(b);
    if (a != b) return error.Mismatch;
    const cells: f64 = @floatFromInt(64 * n);
    std.debug.print("\nthroughput m=64 n={d}: myers {d:.2} GCUPS, derived {d:.2} GCUPS (ops: myers {d}, derived {d})\n", .{
        n,
        cells / @as(f64, @floatFromInt(t_myers)),
        cells / @as(f64, @floatFromInt(t_derived)),
        bitdp.reference.myers_ops,
        K.ops,
    });
}
