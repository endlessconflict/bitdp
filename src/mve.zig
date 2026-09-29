//! Minimal viable experiment: derive kernels, print them, verify 10^6 random
//! pairs against the scalar oracle, and time them against Myers' kernel.

const std = @import("std");
const bitdp = @import("bitdp");

fn printPlan(comptime s: bitdp.Scheme) void {
    const plan = comptime bitdp.derive(s);
    std.debug.print("\nscheme match={d} mismatch={d} gap={d} classes={d}: differences {any}, {d} ops/word, {d} carry chain(s), {s}\n", .{
        s.match, s.mismatch, s.gap, plan.cls.n, plan.vals[0..plan.k], plan.cost, plan.chains, @tagName(plan.method),
    });
    if (plan.len <= 40) for (plan.nodes[0..plan.len], 0..) |n, i| {
        std.debug.print("  r{d:<3} = {s:<5} {d} {d}\n", .{ i, @tagName(n.op), n.a, n.b });
    };
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
        const m = rnd.intRangeAtMost(usize, 1, 62);
        const n = rnd.intRangeAtMost(usize, 0, 256);
        for (p[0..m]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
        for (t[0..n]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
        const want = bitdp.reference.scalar(s, p[0..m], t[0..n], &buf);
        const got = K.distance(p[0..m], t[0..n]);
        if (want != got) {
            std.debug.print("MISMATCH m={d} n={d} want={d} got={d}\n", .{ m, n, want, got });
            return error.Mismatch;
        }
    }
    std.debug.print("  verified {d} random pairs (m<=62, n<=256): bit-exact\n", .{pairs});
}

const bitpal = bitdp.schemes.bitpal;
const general = [_]bitdp.Scheme{ bitdp.schemes.ts_tv, bitdp.schemes.blosum62Linear(4) };

fn now(io: std.Io) i96 {
    return std.Io.Timestamp.now(io, .awake).nanoseconds;
}

/// Derived kernel vs the plain scalar DP (our unvectorized oracle), m = 62.
fn vsScalar(io: std.Io, gpa: std.mem.Allocator, comptime s: bitdp.Scheme, name: []const u8) !void {
    const n: usize = 2_000_000;
    const text = try gpa.alloc(u8, n);
    defer gpa.free(text);
    var prng = std.Random.DefaultPrng.init(11);
    for (text) |*x| x.* = s.alphabet[prng.random().uintLessThan(usize, s.alphabet.len)];
    const pat = text[5000..5062];
    var buf: [65]i64 = undefined;
    var t0 = now(io);
    const a = bitdp.reference.scalar(s, pat, text, &buf);
    const t_scalar = now(io) - t0;
    t0 = now(io);
    const b = bitdp.Kernel(s).distance(pat, text);
    const t_derived = now(io) - t0;
    if (a != b) return error.Mismatch;
    const cells: f64 = @floatFromInt(62 * n);
    std.debug.print("{s}: scalar {d:.2} GCUPS, derived {d:.2} GCUPS ({d} ops/word)\n", .{
        name, cells / @as(f64, @floatFromInt(t_scalar)), cells / @as(f64, @floatFromInt(t_derived)), bitdp.Kernel(s).ops,
    });
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    printPlan(.{});
    printPlan(.{ .mismatch = 2 });
    printPlan(.{ .mismatch = 1, .gap = 2 });
    printPlan(.{ .mismatch = 3, .gap = 2 });
    // BitPAl's published weight sets (score M, I, G), as costs (-M, -I, -G).
    inline for (bitpal) |w| printPlan(w);
    inline for (general) |w| printPlan(w);

    try verify(.{}, 1_000_000);
    try verify(.{ .mismatch = 2 }, 200_000);
    try verify(.{ .mismatch = 1, .gap = 2 }, 200_000);
    try verify(.{ .mismatch = 3, .gap = 2 }, 200_000);
    inline for (bitpal) |w| try verify(w, 100_000);
    inline for (general) |w| try verify(w, 100_000);

    // Throughput: m = 64 against a long random text.
    const gpa = init.gpa;
    const n: usize = 50_000_000;
    const text = try gpa.alloc(u8, n);
    defer gpa.free(text);
    var prng = std.Random.DefaultPrng.init(7);
    for (text) |*x| x.* = "ACGT"[prng.random().int(u2)];
    const pat = text[1000..1062];
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
    const cells: f64 = @floatFromInt(62 * n);
    try vsScalar(io, gpa, bitpal[1], "BitPAl (2,-3,-5)");
    try vsScalar(io, gpa, general[1], "BLOSUM62, gap 4");
    std.debug.print("\nthroughput m=62 n={d}: myers {d:.2} GCUPS, derived {d:.2} GCUPS (ops: myers {d}, derived {d})\n", .{
        n,
        cells / @as(f64, @floatFromInt(t_myers)),
        cells / @as(f64, @floatFromInt(t_derived)),
        bitdp.reference.myers_ops,
        K.ops,
    });
}
